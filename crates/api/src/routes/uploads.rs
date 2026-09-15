//! Upload routes: reference (files or NCBI accession), queries, panel.

use crate::error::{ApiError, ApiResult};
use crate::models::FileDto;
use crate::state::SharedState;
use axum::extract::{Multipart, Path, State};
use axum::Json;
use bactiment_engine::fasta;
use std::io::Write;
use uuid::Uuid;

const MAX_UPLOAD_BYTES: usize = 512 * 1024 * 1024;

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
            bactiment_engine::gff::parse_gff_str(&text).map_err(|e| {
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
    delete_role(&state, project_id, "panel").await?;
    let dto = store_upload(&state, project_id, "panel", &name, data).await?;
    Ok(Json(dto))
}

#[derive(serde::Serialize)]
pub struct PanelFromIdsDto {
    pub file: FileDto,
    pub found: Vec<String>,
    /// Genes fetched from NCBI because the reference lacks them.
    pub from_ncbi: Vec<String>,
    pub missing: Vec<String>,
}

/// POST /projects/{id}/panel/from_ids : a CSV/TSV file of gene identifiers.
/// The panel FASTA is generated automatically: genes present in the
/// reference are extracted from it, the rest are looked up on NCBI.
pub async fn upload_panel_ids(
    State(state): State<SharedState>,
    Path(project_id): Path<i64>,
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
    build_panel(state, project_id, &ids_text, &name).await
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
    Json(body): Json<PanelTextBody>,
) -> ApiResult<Json<PanelFromIdsDto>> {
    ensure_project(&state, project_id).await?;
    if body.text.trim().is_empty() {
        return Err(ApiError::BadRequest(
            "The gene list is empty. Enter gene names separated by commas or new lines.".into(),
        ));
    }
    build_panel(state, project_id, &body.text, "gene list").await
}

async fn build_panel(
    state: SharedState,
    project_id: i64,
    ids_text: &str,
    source_name: &str,
) -> ApiResult<Json<PanelFromIdsDto>> {
    let Some((ref_fasta, ref_gff)) = reference_paths(&state, project_id)? else {
        return Err(ApiError::BadRequest(
            "Please add the reference genome (FASTA + GFF) first: the gene panel is built from it."
                .into(),
        ));
    };
    let ids_text = ids_text.to_string();
    let panel = tokio::task::spawn_blocking(move || {
        bactiment_engine::panel::panel_from_ids(&ref_fasta, &ref_gff, &ids_text)
    })
    .await
    .map_err(|e| ApiError::Internal(format!("The panel could not be built. ({e})")))??;

    // Genes the reference does not carry: fetch them from NCBI by name.
    let mut fasta = panel.fasta;
    let mut from_ncbi = Vec::new();
    let mut missing = panel.missing.clone();
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
                .user_agent("bactiment/0.1")
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

    if panel.found.is_empty() && from_ncbi.is_empty() {
        return Err(ApiError::BadRequest(format!(
            "None of the entries in \u{201c}{source_name}\u{201d} match a gene in the reference annotation, and none could be fetched from NCBI. Check the spelling of the gene names."
        )));
    }

    delete_role(&state, project_id, "panel").await?;
    let dto = store_upload(
        &state,
        project_id,
        "panel",
        "genes_of_interest.fasta",
        fasta.into_bytes(),
    )
    .await?;
    Ok(Json(PanelFromIdsDto {
        file: dto,
        found: panel.found,
        from_ncbi,
        missing,
    }))
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

/// DELETE /projects/{id}/files/{file_id}
pub async fn delete_file(
    State(state): State<SharedState>,
    Path((project_id, file_id)): Path<(i64, i64)>,
) -> ApiResult<&'static str> {
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
    Ok("deleted")
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
