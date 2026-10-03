//! Upload routes: reference (files or NCBI accession), queries, panel.

use crate::error::{ApiError, ApiResult};
use crate::models::FileDto;
use crate::state::SharedState;
use axum::extract::{Multipart, Path, Query, State};
use axum::Json;
use std::io::Write;
use straincompass_engine::fasta;
use uuid::Uuid;

const MAX_UPLOAD_BYTES: usize = 512 * 1024 * 1024;

/// `?append=true` on the panel endpoints: add the new genes to the
/// project's current panel instead of replacing it.
#[derive(serde::Deserialize, Default)]
pub struct PanelMode {
    #[serde(default)]
    pub append: bool,
}

/// Replace the project's panel with `fasta`, or with `fasta` merged into
/// the current panel when appending. A gene already in the panel is
/// replaced by its new sequence, so re-adding a gene fixes it rather
/// than duplicating it. Appending keeps the current panel's file name.
async fn store_panel(
    state: &SharedState,
    project_id: i64,
    append: bool,
    display_name: &str,
    fasta_text: String,
) -> ApiResult<FileDto> {
    let current: Option<(String, String)> = if append {
        let conn = state.db.lock().unwrap();
        conn.query_row(
            "SELECT display_name, stored_name FROM files
             WHERE project_id = ?1 AND role = 'panel' ORDER BY created_at DESC LIMIT 1",
            [project_id],
            |r| Ok((r.get(0)?, r.get(1)?)),
        )
        .ok()
    } else {
        None
    };
    let (name, bytes) = match current {
        Some((cur_name, stored)) => {
            let old = std::fs::read_to_string(state.uploads_dir(project_id).join(stored))
                .unwrap_or_default();
            let new_recs = fasta::parse_fasta_str(&fasta_text).map_err(|e| {
                ApiError::BadRequest(format!("\u{201c}{display_name}\u{201d}: {e}"))
            })?;
            let mut merged = String::new();
            if let Ok(old_recs) = fasta::parse_fasta_str(&old) {
                for r in old_recs
                    .iter()
                    // a gene added again replaces all its old variants
                    .filter(|r| {
                        let gene = straincompass_engine::panel::variant_gene(&r.id);
                        new_recs
                            .iter()
                            .all(|n| straincompass_engine::panel::variant_gene(&n.id) != gene)
                    })
                {
                    push_record(&mut merged, r);
                }
            }
            for r in &new_recs {
                push_record(&mut merged, r);
            }
            (cur_name, merged.into_bytes())
        }
        None => (display_name.to_string(), fasta_text.into_bytes()),
    };
    delete_role(state, project_id, "panel").await?;
    store_upload(state, project_id, "panel", &name, bytes).await
}

fn push_record(out: &mut String, r: &fasta::FastaRecord) {
    out.push('>');
    out.push_str(&r.id);
    if !r.desc.is_empty() {
        out.push(' ');
        out.push_str(&r.desc);
    }
    out.push('\n');
    for chunk in r.seq.chunks(60) {
        out.push_str(&String::from_utf8_lossy(chunk));
        out.push('\n');
    }
}

fn safe_filename(name: &str) -> String {
    let cleaned: String = name
        .chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || matches!(c, '.' | '-' | '_') {
                c
            } else {
                '_'
            }
        })
        .collect();
    if cleaned.is_empty() {
        "file".into()
    } else {
        cleaned
    }
}

/// Validate + store one upload. Returns the DB row.
pub async fn store_upload_public(
    state: &SharedState,
    project_id: i64,
    role: &str,
    display_name: &str,
    bytes: Vec<u8>,
) -> ApiResult<FileDto> {
    store_upload(state, project_id, role, display_name, bytes).await
}

async fn store_upload(
    state: &SharedState,
    project_id: i64,
    role: &str,
    display_name: &str,
    bytes: Vec<u8>,
) -> ApiResult<FileDto> {
    if bytes.len() > MAX_UPLOAD_BYTES {
        return Err(ApiError::BadRequest(
            "This file is too large. The limit is 512 MB.".into(),
        ));
    }
    if bytes.is_empty() {
        return Err(ApiError::BadRequest(
            "The uploaded file is empty. Please check the file and try again.".into(),
        ));
    }
    // Content validation with friendly errors.
    let text = String::from_utf8_lossy(&bytes);
    match role {
        "reference_fasta" | "query" | "panel" => {
            fasta::parse_fasta_str(&text).map_err(|e| {
                ApiError::BadRequest(format!("\u{201c}{display_name}\u{201d}: {}", e))
            })?;
        }
        "reference_gff" => {
            straincompass_engine::gff::parse_gff_str(&text).map_err(|e| {
                ApiError::BadRequest(format!("\u{201c}{display_name}\u{201d}: {}", e))
            })?;
        }
        _ => {}
    }
    let uploads = state.uploads_dir(project_id);
    std::fs::create_dir_all(&uploads)?;
    let stored_name = format!(
        "{}_{}",
        Uuid::new_v4().simple(),
        safe_filename(display_name)
    );
    let path = uploads.join(&stored_name);
    let mut f = std::fs::File::create(&path)?;
    f.write_all(&bytes)?;
    drop(f);
    let created_at = crate::db::now_rfc3339();
    let id = {
        let conn = state.db.lock().unwrap();
        conn.execute(
            "INSERT INTO files (project_id, role, display_name, stored_name, size, created_at)
             VALUES (?1, ?2, ?3, ?4, ?5, ?6)",
            rusqlite::params![
                project_id,
                role,
                display_name,
                stored_name,
                bytes.len() as i64,
                created_at
            ],
        )?;
        conn.last_insert_rowid()
    };
    Ok(FileDto {
        id,
        role: role.into(),
        display_name: display_name.to_string(),
        size: bytes.len() as u64,
        created_at,
    })
}

/// POST /projects/{id}/reference : multipart fields `fasta`, `gff`.
pub async fn upload_reference(
    State(state): State<SharedState>,
    Path(project_id): Path<i64>,
    mut multipart: Multipart,
) -> ApiResult<Json<Vec<FileDto>>> {
    ensure_project(&state, project_id).await?;
    let mut fasta_bytes: Option<(String, Vec<u8>)> = None;
    let mut gff_bytes: Option<(String, Vec<u8>)> = None;
    while let Some(field) = multipart.next_field().await.map_err(|_| {
        ApiError::BadRequest("The upload could not be read. Please try again.".into())
    })? {
        let name = field.name().unwrap_or("").to_string();
        let filename = field.file_name().unwrap_or("file").to_string();
        let data = field.bytes().await.map_err(|_| {
            ApiError::BadRequest("The upload was interrupted. Please try again.".into())
        })?;
        match name.as_str() {
            "fasta" => fasta_bytes = Some((filename, data.to_vec())),
            "gff" => gff_bytes = Some((filename, data.to_vec())),
            _ => {}
        }
    }
    let (fname, fdata) = fasta_bytes.ok_or_else(|| {
        ApiError::BadRequest(
            "Please provide both the reference genome (FASTA) and its annotation (GFF).".into(),
        )
    })?;
    let (gname, gdata) = gff_bytes.ok_or_else(|| {
        ApiError::BadRequest(
            "Please provide both the reference genome (FASTA) and its annotation (GFF).".into(),
        )
    })?;
    let mut created = Vec::new();
    {
        // replace any previous reference
        delete_role(&state, project_id, "reference_fasta").await?;
        delete_role(&state, project_id, "reference_gff").await?;
    }
    created.push(store_upload(&state, project_id, "reference_fasta", &fname, fdata).await?);
    created.push(store_upload(&state, project_id, "reference_gff", &gname, gdata).await?);
    Ok(Json(created))
}

/// POST /projects/{id}/queries : one or many FASTA files.
pub async fn upload_queries(
    State(state): State<SharedState>,
    Path(project_id): Path<i64>,
    mut multipart: Multipart,
) -> ApiResult<Json<Vec<FileDto>>> {
    ensure_project(&state, project_id).await?;
    let mut created = Vec::new();
    while let Some(field) = multipart.next_field().await.map_err(|_| {
        ApiError::BadRequest("The upload could not be read. Please try again.".into())
    })? {
        let filename = field.file_name().unwrap_or("").to_string();
        let data = field.bytes().await.map_err(|_| {
            ApiError::BadRequest("The upload was interrupted. Please try again.".into())
        })?;
        if filename.is_empty() {
            continue;
        }
        created.push(store_upload(&state, project_id, "query", &filename, data.to_vec()).await?);
    }
    if created.is_empty() {
        return Err(ApiError::BadRequest(
            "No query files were received. Please choose at least one FASTA file.".into(),
        ));
    }
    Ok(Json(created))
}

/// POST /projects/{id}/panel : one FASTA file.
pub async fn upload_panel(
    State(state): State<SharedState>,
    Path(project_id): Path<i64>,
    Query(mode): Query<PanelMode>,
    mut multipart: Multipart,
) -> ApiResult<Json<FileDto>> {
    ensure_project(&state, project_id).await?;
    let mut got: Option<(String, Vec<u8>)> = None;
    while let Some(field) = multipart.next_field().await.map_err(|_| {
        ApiError::BadRequest("The upload could not be read. Please try again.".into())
    })? {
        let filename = field.file_name().unwrap_or("").to_string();
        let data = field.bytes().await.map_err(|_| {
            ApiError::BadRequest("The upload was interrupted. Please try again.".into())
        })?;
        if filename.is_empty() {
            continue;
        }
        got = Some((filename, data.to_vec()));
    }
    let (name, data) = got.ok_or_else(|| {
        ApiError::BadRequest("No gene panel file was received. Please choose a FASTA file.".into())
    })?;
    let text = String::from_utf8_lossy(&data).into_owned();
    let dto = store_panel(&state, project_id, mode.append, &name, text).await?;
    Ok(Json(dto))
}

#[derive(serde::Serialize)]
pub struct PanelFromIdsDto {
    pub file: FileDto,
    pub found: Vec<String>,
    /// Genes fetched from NCBI because the reference lacks them.
    pub from_ncbi: Vec<String>,
    /// Genes taken from the curated catalogs (AMRFinderPlus, VFDB), with
    /// the entry and product they were taken from.
    #[serde(default)]
    pub from_catalog: Vec<String>,
    /// Genes taken from the local reference library, with their variants.
    #[serde(default)]
    pub from_library: Vec<String>,
    /// The library used ("Listeria library 2026-10-02"), if any.
    #[serde(default)]
    pub library: Option<String>,
    /// Genes taken from the reference for which the curated catalogs hold
    /// a different, organism-specific gene of the same name.
    #[serde(default)]
    /// The notes on the panel, one group per gene they concern.
    pub hints: Vec<GeneNotes>,
    pub missing: Vec<String>,
}

/// POST /projects/{id}/panel/from_ids : a CSV/TSV file of gene identifiers.
/// The panel FASTA is generated automatically: genes present in the
/// reference are extracted from it, the rest are looked up on NCBI.
pub async fn upload_panel_ids(
    State(state): State<SharedState>,
    Path(project_id): Path<i64>,
    Query(mode): Query<PanelMode>,
    mut multipart: Multipart,
) -> ApiResult<Json<PanelFromIdsDto>> {
    ensure_project(&state, project_id).await?;
    let mut got: Option<(String, Vec<u8>)> = None;
    while let Some(field) = multipart.next_field().await.map_err(|_| {
        ApiError::BadRequest("The upload could not be read. Please try again.".into())
    })? {
        let filename = field.file_name().unwrap_or("genes.txt").to_string();
        let data = field.bytes().await.map_err(|_| {
            ApiError::BadRequest("The upload was interrupted. Please try again.".into())
        })?;
        if filename.is_empty() {
            continue;
        }
        got = Some((filename, data.to_vec()));
    }
    let (name, data) = got.ok_or_else(|| {
        ApiError::BadRequest(
            "No gene list file was received. Please choose a CSV or TSV file.".into(),
        )
    })?;
    let ids_text = String::from_utf8_lossy(&data).to_string();
    if ids_text.lines().all(|l| l.trim().is_empty()) {
        return Err(ApiError::BadRequest(
            "This file appears to be empty. Please check the file and try again.".into(),
        ));
    }
    build_panel(state, project_id, &ids_text, &name, mode.append).await
}

#[derive(serde::Deserialize)]
pub struct PanelTextBody {
    pub text: String,
}

/// POST /projects/{id}/panel/from_text : gene identifiers pasted as text
/// (commas, spaces or new lines as separators).
pub async fn upload_panel_text(
    State(state): State<SharedState>,
    Path(project_id): Path<i64>,
    Query(mode): Query<PanelMode>,
    Json(body): Json<PanelTextBody>,
) -> ApiResult<Json<PanelFromIdsDto>> {
    ensure_project(&state, project_id).await?;
    if body.text.trim().is_empty() {
        return Err(ApiError::BadRequest(
            "The gene list is empty. Enter gene names separated by commas or new lines.".into(),
        ));
    }
    build_panel(state, project_id, &body.text, "gene list", mode.append).await
}

async fn build_panel(
    state: SharedState,
    project_id: i64,
    ids_text: &str,
    source_name: &str,
    append: bool,
) -> ApiResult<Json<PanelFromIdsDto>> {
    let Some((ref_fasta, ref_gff)) = reference_paths(&state, project_id)? else {
        return Err(ApiError::BadRequest(
            "Please add the reference genome (FASTA + GFF) first: the gene panel is built from it."
                .into(),
        ));
    };
    let ids_text = ids_text.to_string();
    let chooser_lib = crate::routes::library::for_project(&state, project_id);
    let chooser_organism: String = {
        let conn = state.db.lock().unwrap();
        conn.query_row(
            "SELECT organism FROM projects WHERE id = ?1",
            [project_id],
            |r| r.get::<_, Option<String>>(0),
        )
        .ok()
        .flatten()
        .unwrap_or_default()
    };
    let chooser_vfdb = state.data_dir.join("db").join("VFDB_setA_nt.fas");
    let chooser_work = state.cache_dir(project_id).join("library");
    let panel = tokio::task::spawn_blocking(move || {
        let catalogs = load_catalogs(chooser_vfdb);
        let blastn = straincompass_engine::tools::find_optional("blastn");
        let choose = |name: &str, cands: &[straincompass_engine::panel::Candidate]| {
            choose_gene(
                name,
                cands,
                chooser_lib.as_ref(),
                blastn.as_deref(),
                &catalogs,
                &chooser_organism,
                &chooser_work,
            )
        };
        straincompass_engine::panel::panel_from_ids_with(&ref_fasta, &ref_gff, &ids_text, &choose)
    })
    .await
    .map_err(|e| ApiError::Internal(format!("The panel could not be built. ({e})")))??;

    // Genes the reference does not carry: first the local reference
    // library of the project's genus (every variant of the name, from the
    // genus' own genomes), then the curated catalogs the resistance/
    // virulence screen uses (sequence-curated, organism-aware names), then
    // NCBI by name.
    let mut fasta = panel.fasta;
    let mut from_ncbi = Vec::new();
    let mut missing = panel.missing.clone();
    let mut from_catalog = Vec::new();
    let mut from_library = Vec::new();
    let mut library_hints = Vec::new();
    let mut library_title = None;
    if let Some(lib) = crate::routes::library::for_project(&state, project_id) {
        library_title = Some(lib.title());
        let wanted = missing.clone();
        let found = panel.found.clone();
        let ref_panel = fasta.clone();
        let ref_genome = reference_paths(&state, project_id)?.map(|(f, _)| f);
        let work = state.cache_dir(project_id).join("library");
        let looked = tokio::task::spawn_blocking(move || {
            let blastn = straincompass_engine::tools::find_optional("blastn");
            let blastx = straincompass_engine::tools::find_optional("blastx");
            crate::routes::library::panel_lookup(
                &lib,
                blastn.as_deref(),
                blastx.as_deref(),
                &wanted,
                &found,
                &ref_panel,
                ref_genome.as_deref(),
                &work,
            )
        })
        .await
        .map_err(|e| {
            ApiError::Internal(format!("The reference library could not be read. ({e})"))
        })?;
        fasta.push_str(&looked.records);
        from_library = looked.notes;
        library_hints = looked.hints;
        library_hints.splice(0..0, panel.ambiguous.iter().cloned());
        missing = looked.unresolved;
    }
    let organism: String = {
        let conn = state.db.lock().unwrap();
        conn.query_row(
            "SELECT organism FROM projects WHERE id = ?1",
            [project_id],
            |r| r.get::<_, Option<String>>(0),
        )
        .ok()
        .flatten()
        .unwrap_or_default()
    };
    let vfdb_default = state.data_dir.join("db").join("VFDB_setA_nt.fas");
    let wanted = missing.clone();
    let ref_panel = fasta.clone();
    let found = panel.found.clone();
    let organism2 = organism.clone();
    let (resolved, unresolved, mut hints) = tokio::task::spawn_blocking(move || {
        catalog_lookup(&wanted, &found, &ref_panel, &organism2, vfdb_default)
    })
    .await
    .map_err(|e| ApiError::Internal(format!("The gene catalogs could not be read. ({e})")))?;
    for (record, note) in resolved {
        fasta.push_str(&record);
        if !note.is_empty() {
            from_catalog.push(note);
        }
    }
    missing = unresolved;
    if library_title.is_none() {
        // without a library the ambiguity notes still need saying
        library_hints.splice(0..0, panel.ambiguous.iter().cloned());
    }
    hints.splice(0..0, library_hints);
    if !missing.is_empty() {
        let (organism, api_key) = {
            let conn = state.db.lock().unwrap();
            let organism: String = conn
                .query_row(
                    "SELECT organism FROM projects WHERE id = ?1",
                    [project_id],
                    |r| r.get(0),
                )
                .unwrap_or_default();
            let key = crate::db::get_setting(&conn, "ncbi_api_key")?;
            (organism, key)
        };
        let genus = organism.split_whitespace().next().unwrap_or("").to_string();
        if !genus.is_empty() {
            let c = reqwest::Client::builder()
                .timeout(std::time::Duration::from_secs(60))
                .connect_timeout(std::time::Duration::from_secs(20))
                .user_agent("straincompass/0.1")
                .build()
                .map_err(|e| ApiError::Internal(format!("NCBI could not be reached. ({e})")))?;
            let mut unresolved: Vec<String> = Vec::new();
            for line in &missing {
                let mut got_one = false;
                // a GenBank accession pinned in the list itself:
                // "qacH (HF565366.1)" or "emrC (CP038643.1:1496-1882 rev)"
                if let Some(spec) = crate::routes::ncbi::parse_accession_spec(line) {
                    if let Some(g) =
                        crate::routes::ncbi::fetch_gene_by_accession(&c, api_key.as_deref(), &spec)
                            .await
                    {
                        fasta.push_str(&g.record);
                        from_ncbi.push(format!(
                            "{} ({})",
                            spec.name.clone().unwrap_or_else(|| spec.accession.clone()),
                            g.source
                        ));
                        got_one = true;
                    }
                }
                if got_one {
                    continue;
                }
                // "pva (lmo0446)": try the parenthesized locus tag first
                // (a specific genome), then the bare symbol as fallback
                let mut tokens: Vec<String> = Vec::new();
                for part in line.split([',', '\t', ';']) {
                    for tok in part.split_whitespace() {
                        let bare = tok
                            .trim_start_matches(['(', '['])
                            .trim_end_matches([')', ']']);
                        if bare.is_empty() || bare.starts_with('#') {
                            continue;
                        }
                        if tok.starts_with('(') || tok.starts_with('[') {
                            tokens.insert(0, bare.to_string());
                        } else {
                            tokens.push(bare.to_string());
                        }
                    }
                }
                for name in tokens {
                    if let Some(g) = crate::routes::ncbi::fetch_gene(
                        &c,
                        api_key.as_deref(),
                        &organism,
                        &genus,
                        &name,
                    )
                    .await
                    {
                        fasta.push_str(&g.record);
                        from_ncbi.push(format!("{} ({})", name, g.source));
                        got_one = true;
                        break;
                    }
                }
                if !got_one {
                    unresolved.push(line.clone());
                }
            }
            missing = unresolved;
        }
    }

    if panel.found.is_empty()
        && from_ncbi.is_empty()
        && from_catalog.is_empty()
        && from_library.is_empty()
    {
        let hints: Vec<String> = missing.iter().filter_map(|l| accession_hint(l)).collect();
        let mut msg = format!(
            "None of the entries in \u{201c}{source_name}\u{201d} match a gene in the reference annotation, and none could be fetched from NCBI. Check the spelling of the gene names."
        );
        for h in hints {
            msg.push(' ');
            msg.push_str(&h);
        }
        return Err(ApiError::BadRequest(msg));
    }

    // the sources each numbered their own variants: number them again,
    // dropping a sequence a gene already holds
    if let Ok(recs) = fasta::parse_fasta_str(&fasta) {
        fasta = straincompass_engine::panel::renumber_variants(&recs);
    }

    // variant sets: operon partners from the records the panel's genes
    // came from, and a word on sequences from another genus
    let api_key = {
        let conn = state.db.lock().unwrap();
        crate::db::get_setting(&conn, "ncbi_api_key")?
    };
    let genus = organism.split_whitespace().next().unwrap_or("");
    let lib = crate::routes::library::for_project(&state, project_id);
    let (extra, notes, summary) =
        crate::routes::ncbi::panel_variant_sets(&fasta, genus, api_key.as_deref(), lib.as_ref())
            .await;
    fasta.push_str(&extra);
    hints.extend(notes);
    // the plain list of a gene's variants, only where nothing above
    // already tells them
    let told = |g: &str| {
        let lead = format!("{g}: ");
        hints.iter().any(|h| h.starts_with(&lead))
            || from_library
                .iter()
                .any(|n| n.starts_with(&format!("{g} \u{2190}")))
    };
    let summary: Vec<String> = summary
        .into_iter()
        .filter(|(g, _)| !told(g))
        .map(|(_, line)| line)
        .collect();
    hints.extend(summary);
    let hints = by_gene(hints);

    let dto = store_panel(&state, project_id, append, "genes_of_interest.fasta", fasta).await?;
    Ok(Json(PanelFromIdsDto {
        file: dto,
        found: panel.found,
        from_ncbi,
        from_catalog,
        from_library,
        library: library_title,
        hints,
        missing,
    }))
}

/// A gene's notes on the panel build, told together.
#[derive(Debug, serde::Serialize, PartialEq)]
pub struct GeneNotes {
    /// Empty for notes on no one gene.
    pub gene: String,
    pub notes: Vec<String>,
}

/// Notes grouped by the gene each leads with ("prfA: Taken ..."), in the
/// order the genes first come up, the lead taken off.
fn by_gene(hints: Vec<String>) -> Vec<GeneNotes> {
    let mut out: Vec<GeneNotes> = Vec::new();
    for h in hints {
        let (gene, text) = match h.split_once(": ") {
            Some((g, t)) if !g.is_empty() && !g.contains(char::is_whitespace) => {
                (g.to_string(), t.to_string())
            }
            _ => (String::new(), h),
        };
        let mut c = text.chars();
        let text = match c.next() {
            Some(f) => f.to_uppercase().chain(c).collect(),
            None => continue,
        };
        match out.iter_mut().find(|g| g.gene == gene) {
            Some(g) => {
                if !g.notes.contains(&text) {
                    g.notes.push(text)
                }
            }
            None => out.push(GeneNotes {
                gene,
                notes: vec![text],
            }),
        }
    }
    out
}

/// GET /projects/{id}/files
pub async fn list_files(
    State(state): State<SharedState>,
    Path(project_id): Path<i64>,
) -> ApiResult<Json<Vec<FileDto>>> {
    ensure_project(&state, project_id).await?;
    let conn = state.db.lock().unwrap();
    let mut stmt = conn.prepare(
        "SELECT id, role, display_name, size, created_at FROM files
         WHERE project_id = ?1 ORDER BY role, created_at",
    )?;
    let rows = stmt
        .query_map([project_id], |r| {
            Ok(FileDto {
                id: r.get(0)?,
                role: r.get(1)?,
                display_name: r.get(2)?,
                size: r.get::<_, i64>(3)? as u64,
                created_at: r.get(4)?,
            })
        })?
        .collect::<std::result::Result<Vec<_>, _>>()?;
    Ok(Json(rows))
}

/// GET /projects/{id}/files/{file_id} : the file exactly as it was
/// uploaded (or built, for a panel), under its original name.
pub async fn download_file(
    State(state): State<SharedState>,
    Path((project_id, file_id)): Path<(i64, i64)>,
) -> ApiResult<axum::response::Response> {
    use axum::response::IntoResponse;
    ensure_project(&state, project_id).await?;
    let row: Option<(String, String)> = {
        let conn = state.db.lock().unwrap();
        conn.query_row(
            "SELECT display_name, stored_name FROM files WHERE id = ?1 AND project_id = ?2",
            [file_id, project_id],
            |r| Ok((r.get(0)?, r.get(1)?)),
        )
        .ok()
    };
    let Some((display_name, stored_name)) = row else {
        return Err(ApiError::NotFound(
            "This file does not exist (anymore).".into(),
        ));
    };
    let bytes = tokio::fs::read(state.uploads_dir(project_id).join(stored_name))
        .await
        .map_err(|_| ApiError::NotFound("This file is missing on the server.".into()))?;
    let body = axum::body::Body::from(bytes);
    let disposition = format!("attachment; filename=\"{}\"", safe_filename(&display_name));
    Ok((
        [
            (
                axum::http::header::CONTENT_TYPE,
                "application/octet-stream".to_string(),
            ),
            (axum::http::header::CONTENT_DISPOSITION, disposition),
        ],
        body,
    )
        .into_response())
}

/// DELETE /projects/{id}/files/{file_id}
pub async fn delete_file(
    State(state): State<SharedState>,
    Path((project_id, file_id)): Path<(i64, i64)>,
) -> ApiResult<Json<serde_json::Value>> {
    ensure_project(&state, project_id).await?;
    let stored_name: Option<String> = {
        let conn = state.db.lock().unwrap();
        conn.query_row(
            "SELECT stored_name FROM files WHERE id = ?1 AND project_id = ?2",
            [file_id, project_id],
            |r| r.get(0),
        )
        .ok()
    };
    let Some(stored_name) = stored_name else {
        return Err(ApiError::NotFound(
            "This file does not exist (anymore).".into(),
        ));
    };
    {
        let conn = state.db.lock().unwrap();
        conn.execute("DELETE FROM files WHERE id = ?1", [file_id])?;
    }
    let path = state.uploads_dir(project_id).join(stored_name);
    let _ = std::fs::remove_file(path);
    Ok(Json(serde_json::json!({ "status": "deleted" })))
}

pub async fn delete_role(state: &SharedState, project_id: i64, role: &str) -> ApiResult<()> {
    let names: Vec<String> = {
        let conn = state.db.lock().unwrap();
        let mut stmt =
            conn.prepare("SELECT stored_name FROM files WHERE project_id = ?1 AND role = ?2")?;
        let rows = stmt
            .query_map(rusqlite::params![project_id, role], |r| r.get(0))?
            .collect::<std::result::Result<Vec<_>, _>>()?;
        rows
    };
    if !names.is_empty() {
        let conn = state.db.lock().unwrap();
        conn.execute(
            "DELETE FROM files WHERE project_id = ?1 AND role = ?2",
            rusqlite::params![project_id, role],
        )?;
        drop(conn);
        for n in names {
            let _ = std::fs::remove_file(state.uploads_dir(project_id).join(n));
        }
    }
    Ok(())
}

pub async fn ensure_project(state: &SharedState, project_id: i64) -> ApiResult<()> {
    let exists = {
        let conn = state.db.lock().unwrap();
        conn.query_row("SELECT 1 FROM projects WHERE id = ?1", [project_id], |_| {
            Ok(())
        })
        .is_ok()
    };
    if !exists {
        return Err(ApiError::NotFound(
            "This project does not exist (anymore).".into(),
        ));
    }
    Ok(())
}

/// Resolve the current reference (fasta, gff) of a project.
pub fn reference_paths(
    state: &SharedState,
    project_id: i64,
) -> ApiResult<Option<(std::path::PathBuf, std::path::PathBuf)>> {
    let conn = state.db.lock().unwrap();
    let mut stmt = conn.prepare(
        "SELECT role, stored_name FROM files WHERE project_id = ?1 AND role IN ('reference_fasta','reference_gff')",
    )?;
    let rows = stmt
        .query_map([project_id], |r| {
            Ok((r.get::<_, String>(0)?, r.get::<_, String>(1)?))
        })?
        .collect::<std::result::Result<Vec<_>, _>>()?;
    drop(stmt);
    let mut fasta = None;
    let mut gff = None;
    for (role, stored) in rows {
        match role.as_str() {
            "reference_fasta" => fasta = Some(state.uploads_dir(project_id).join(stored)),
            "reference_gff" => gff = Some(state.uploads_dir(project_id).join(stored)),
            _ => {}
        }
    }
    Ok(match (fasta, gff) {
        (Some(f), Some(g)) => Some((f, g)),
        _ => None,
    })
}

/// Why a "name (ACC:start-end)" entry could not be used, when the part in
/// brackets is not a GenBank accession (a copied placeholder, a typo).
fn accession_hint(line: &str) -> Option<String> {
    let open = line.find('(')?;
    let close = line.rfind(')')?;
    let inner = line.get(open + 1..close)?;
    let acc = inner.split(':').next()?.trim();
    if acc.is_empty() || crate::routes::ncbi::parse_accession_spec(line).is_some() {
        return None;
    }
    Some(format!(
        "In \u{201c}{}\u{201d}, \u{201c}{acc}\u{201d} is not a GenBank accession: write the record's real accession, e.g. emrC (CP038643.1:1496-1882 rev).",
        line.trim()
    ))
}

/// Resolve plain gene names (no pinned accession) from the curated
/// catalogs: (FASTA record, note) for each one found, and the entries
/// still unresolved.
/// The curated catalogs the screen uses (AMRFinderPlus CDS, VFDB), as far
/// as this server has them.
fn load_catalogs(vfdb_default: std::path::PathBuf) -> Vec<straincompass_engine::catalog::Catalog> {
    use straincompass_engine::catalog::Catalog;
    let cfg = straincompass_engine::screen::ScreenConfig::discover(Some(vfdb_default));
    let mut catalogs = Vec::new();
    if let Some(c) = cfg
        .amr_db
        .as_ref()
        .and_then(|d| Catalog::from_file(&d.join("AMR_CDS.fa"), false))
    {
        catalogs.push(c);
    }
    if let Some(c) = cfg.vfdb.as_ref().and_then(|p| Catalog::from_file(p, true)) {
        catalogs.push(c);
    }
    catalogs
}

/// Which of the reference's genes of one name the user means: the one a
/// curated entry of that name is. The library knows which variant group
/// each entry matches (by protein); without it, the catalogs' own DNA is
/// compared. None leaves the first in genome order.
fn choose_gene(
    name: &str,
    cands: &[straincompass_engine::panel::Candidate],
    lib: Option<&straincompass_engine::library::Library>,
    blastn: Option<&std::path::Path>,
    catalogs: &[straincompass_engine::catalog::Catalog],
    organism: &str,
    work: &std::path::Path,
) -> Option<(usize, String)> {
    if let (Some(lib), Some(blastn)) = (lib, blastn) {
        let curated = lib.curated_groups(name).unwrap_or_default();
        if !curated.is_empty() {
            let seqs: Vec<(String, Vec<u8>)> = cands
                .iter()
                .enumerate()
                .map(|(i, c)| (format!("c{i}"), c.seq.clone()))
                .collect();
            if let Ok(groups) = lib.groups_of(blastn, &seqs, work) {
                for (group, entry) in &curated {
                    if let Some(i) = (0..cands.len()).find(|i| {
                        groups
                            .get(&format!("c{i}"))
                            .is_some_and(|g| g.contains(group))
                    }) {
                        return Some((i, format!("it is what {entry} means")));
                    }
                }
            }
        }
    }
    // the catalogs' own sequences: the candidate sharing the most 21-mers
    // with an entry of that name, if it shares a fair part (another
    // strain's allele shares most; an unrelated gene none)
    let organism = (!organism.trim().is_empty()).then_some(organism);
    let mut best: Option<(f64, usize, String)> = None;
    for g in straincompass_engine::catalog::lookup_all(name, organism, catalogs) {
        for (i, c) in cands.iter().enumerate() {
            let share = straincompass_engine::panel_variants::kmer_share(&c.seq, &g.seq);
            if share >= 0.3 && best.as_ref().is_none_or(|b| share > b.0) {
                best = Some((
                    share,
                    i,
                    format!("it is what {} {} means", g.source, g.symbol),
                ));
            }
        }
    }
    best.map(|(_, i, why)| (i, why))
}

fn catalog_lookup(
    wanted: &[String],
    found: &[String],
    ref_panel: &str,
    organism: &str,
    vfdb_default: std::path::PathBuf,
) -> (Vec<(String, String)>, Vec<String>, Vec<String>) {
    use straincompass_engine::catalog::lookup_all;
    use straincompass_engine::panel::variant_id;
    let catalogs = load_catalogs(vfdb_default);
    let organism = (!organism.trim().is_empty()).then_some(organism);
    let mut resolved = Vec::new();
    let mut unresolved = Vec::new();
    for line in wanted {
        // a pinned accession is the user's explicit choice: leave it to NCBI
        let plain = crate::routes::ncbi::parse_accession_spec(line).is_none();
        let name = line
            .split_whitespace()
            .next()
            .unwrap_or("")
            .trim_matches(|c| c == '(' || c == ')' || c == '[' || c == ']');
        let all = if plain {
            lookup_all(name, organism, &catalogs)
        } else {
            Vec::new()
        };
        match all.first().cloned() {
            Some(g) => {
                // every sequence filed under the entry, as variants
                let mut rec = String::new();
                for (i, v) in all.iter().enumerate() {
                    rec.push_str(&catalog_record(&variant_id(name, i + 1), v));
                }
                let note = format!(
                    "{name} \u{2190} {} {}: {}{}",
                    g.source,
                    g.symbol,
                    g.product,
                    if g.origin.is_empty() {
                        String::new()
                    } else {
                        format!(" ({})", g.origin)
                    }
                );
                resolved.push((rec, note));
            }
            None => unresolved.push(line.clone()),
        }
    }
    // reference genes the catalogs know in another version: searched
    // beside the reference's own copy, which stays the first variant
    let mut hints = Vec::new();
    if let Ok(recs) = straincompass_engine::fasta::parse_fasta_str(ref_panel) {
        for name in found {
            let Some(rec) = recs.iter().find(|r| &r.id == name) else {
                continue;
            };
            if let Some(v) =
                straincompass_engine::catalog::organism_variant(name, organism, &rec.seq, &catalogs)
            {
                let mut extra = String::new();
                let others = lookup_all(name, organism, &catalogs);
                let others: Vec<_> = others
                    .iter()
                    .filter(|g| !g.seq.eq_ignore_ascii_case(&rec.seq))
                    .collect();
                for (i, g) in others.iter().enumerate() {
                    extra.push_str(&catalog_record(&variant_id(name, i + 2), g));
                }
                resolved.push((extra, String::new()));
                hints.push(format!(
                    "{name}: Also searched, beside your reference's copy: the different sequence {} files under {} ({}).",
                    v.source, v.symbol, v.product
                ));
            }
        }
    }
    (resolved, unresolved, hints)
}

/// A catalog gene as a panel record, its header note naming where it is
/// from: "cadA__v2 AMRFinderPlus cadA_Lm: product [AP022822.1:2608314-2610431]".
fn catalog_record(id: &str, g: &straincompass_engine::catalog::CatalogGene) -> String {
    let mut rec = format!(
        ">{id} {} {}: {} [{}]\n",
        g.source, g.symbol, g.product, g.origin
    );
    for chunk in g.seq.chunks(60) {
        rec.push_str(&String::from_utf8_lossy(chunk));
        rec.push('\n');
    }
    rec
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn notes_are_told_once_per_gene() {
        let g = by_gene(vec![
            "prfA: Taken from your reference: x.".into(),
            "hly: Also searched, beside your reference's copy: lso.".into(),
            "prfA: not searched, sharing only the name: RF1.".into(),
            "prfA: not searched, sharing only the name: RF1.".into(),
            "The sequence from AB1.1 is not from Listeria.".into(),
        ]);
        let notes = |gene: &str, notes: &[&str]| GeneNotes {
            gene: gene.into(),
            notes: notes.iter().map(|n| n.to_string()).collect(),
        };
        assert_eq!(
            g,
            [
                notes(
                    "prfA",
                    &[
                        "Taken from your reference: x.",
                        "Not searched, sharing only the name: RF1."
                    ]
                ),
                notes(
                    "hly",
                    &["Also searched, beside your reference's copy: lso."]
                ),
                notes("", &["The sequence from AB1.1 is not from Listeria."]),
            ]
        );
    }
}
