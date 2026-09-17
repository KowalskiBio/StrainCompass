//! The job runner: executes analysis runs in the background.

use crate::error::ApiResult;
use crate::state::SharedState;
use std::io::Write;
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;
use straincompass_engine::pipeline::{
    self, ComparisonInputs, ComparisonResult, QueryAlignmentSource, WorkDirs,
};
use straincompass_engine::tools::ToolPaths;
use straincompass_types::{MatrixRow, RunParams};

pub fn spawn_run(state: SharedState, run_id: i64) {
    tokio::spawn(async move {
        let result = execute_run(&state, run_id).await;
        if let Err(e) = result {
            tracing::error!("run {run_id} failed: {e}");
            let conn = state.db.lock().unwrap();
            let _ = conn.execute(
                "UPDATE runs SET status = 'failed', error = ?2, step = NULL,
                 finished_at = strftime('%Y-%m-%dT%H:%M:%fZ','now') WHERE id = ?1",
                rusqlite::params![run_id, e.to_string()],
            );
            let _ = conn.execute(
                "INSERT INTO run_logs (run_id, seq, line) VALUES (?1, -1, ?2)",
                rusqlite::params![run_id, format!("Run failed: {}", e)],
            );
            drop(conn);
        }
    });
}

#[derive(Clone)]
struct RunCtx {
    state: SharedState,
    run_id: i64,
}

impl RunCtx {
    fn log(&self, line: &str) {
        if let Ok(conn) = self.state.db.lock() {
            let seq: i64 = conn
                .query_row(
                    "SELECT COALESCE(MAX(seq), -1) + 1 FROM run_logs WHERE run_id = ?1",
                    [self.run_id],
                    |r| r.get(0),
                )
                .unwrap_or(0);
            let _ = conn.execute(
                "INSERT OR REPLACE INTO run_logs (run_id, seq, line) VALUES (?1, ?2, ?3)",
                rusqlite::params![self.run_id, seq, line],
            );
        }
        let dir = self.state.data_dir.join("runlogs");
        let _ = std::fs::create_dir_all(&dir);
        let _ = std::fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(dir.join(format!("{}.log", self.run_id)))
            .and_then(|mut f| writeln!(f, "{line}"));
    }
    fn set_step(&self, step: &str) {
        if let Ok(conn) = self.state.db.lock() {
            let _ = conn.execute(
                "UPDATE runs SET step = ?2 WHERE id = ?1",
                rusqlite::params![self.run_id, step],
            );
        }
    }
}

async fn execute_run(state: &SharedState, run_id: i64) -> ApiResult<()> {
    let ctx = RunCtx {
        state: state.clone(),
        run_id,
    };
    // 1. load the run
    let (project_id, params_json, query_ids): (i64, String, String) = {
        let conn = state.db.lock().unwrap();
        conn.query_row(
            "SELECT project_id, params_json, query_ids FROM runs WHERE id = ?1",
            [run_id],
            |r| Ok((r.get(0)?, r.get(1)?, r.get(2)?)),
        )
        .map_err(|_| crate::error::ApiError::NotFound("This run does not exist.".into()))?
    };
    let params: RunParams = serde_json::from_str(&params_json).map_err(|_| {
        crate::error::ApiError::Internal("Stored parameters are unreadable.".into())
    })?;
    let query_file_ids: Vec<i64> = serde_json::from_str(&query_ids)
        .map_err(|_| crate::error::ApiError::Internal("Stored query list is unreadable.".into()))?;

    {
        let conn = state.db.lock().unwrap();
        conn.execute(
            "UPDATE runs SET status = 'running', started_at = strftime('%Y-%m-%dT%H:%M:%fZ','now'),
             step = 'Preparing the comparison' WHERE id = ?1",
            [run_id],
        )?;
    }
    ctx.log("Preparing the comparison.");

    // 2. resolve inputs
    let (ref_fasta, ref_gff) = crate::routes::uploads::reference_paths(state, project_id)?
        .ok_or_else(|| {
            crate::error::ApiError::BadRequest(
                "Please add a reference genome (FASTA + GFF) before comparing.".into(),
            )
        })?;
    let panel_path: Option<(String, PathBuf)> = {
        let conn = state.db.lock().unwrap();
        conn.query_row(
            "SELECT display_name, stored_name FROM files WHERE project_id = ?1 AND role = 'panel'",
            [project_id],
            |r| {
                let display: String = r.get(0)?;
                let stored: String = r.get(1)?;
                Ok((display, state.uploads_dir(project_id).join(stored)))
            },
        )
        .ok()
    };
    let mut queries: Vec<(i64, String, PathBuf)> = Vec::new();
    {
        let conn = state.db.lock().unwrap();
        for qid in &query_file_ids {
            let row: Option<(String, String)> = conn
                .query_row(
                    "SELECT display_name, stored_name FROM files WHERE id = ?1 AND role = 'query'",
                    [qid],
                    |r| Ok((r.get(0)?, r.get(1)?)),
                )
                .ok();
            if let Some((name, stored)) = row {
                queries.push((*qid, name, state.uploads_dir(project_id).join(stored)));
            }
        }
    }
    if queries.is_empty() {
        return Err(crate::error::ApiError::BadRequest(
            "None of the selected query genomes exist anymore.".into(),
        ));
    }

    let _tools = ToolPaths::discover().map_err(|e| {
        crate::error::ApiError::Internal(format!(
            "The analysis tools are not installed on the server. {}",
            e
        ))
    })?;

    // 3. stage the reference (sanitized), cached by content hash
    let ref_dir = state.project_dir(project_id).join("reference");
    std::fs::create_dir_all(&ref_dir)?;
    let staged_fa = ref_dir.join("ref.fa");
    let staged_gff = ref_dir.join("ref.gff");
    let ref_hash = pipeline::file_hash(&ref_fasta)?;
    let gff_hash = pipeline::file_hash(&ref_gff)?;
    let marker = ref_dir.join("staged.txt");
    let staged_ok = std::fs::read_to_string(&marker)
        .map(|m| m.trim() == format!("{ref_hash} {gff_hash}"))
        .unwrap_or(false);
    if !staged_ok {
        ctx.set_step("Reading the reference genome");
        ctx.log("Reading and preparing the reference genome.");
        pipeline::sanitize_into(&ref_fasta, &staged_fa)?;
        std::fs::copy(&ref_gff, &staged_gff)?;
        std::fs::write(&marker, format!("{ref_hash} {gff_hash}"))?;
        // check that GFF seqids exist in the fasta
        let genes = straincompass_engine::gff::parse_gff(&staged_gff)?;
        let recs = straincompass_engine::fasta::parse_fasta(&staged_fa)?;
        let ids: std::collections::HashSet<String> = recs.iter().map(|r| r.id.clone()).collect();
        let missing: Vec<String> = genes
            .iter()
            .filter(|g| !ids.contains(&g.seqid))
            .map(|g| g.seqid.clone())
            .collect::<std::collections::HashSet<_>>()
            .into_iter()
            .collect();
        if !missing.is_empty() {
            return Err(crate::error::ApiError::BadRequest(format!(
                "The annotation (GFF) mentions sequences that are not in the reference FASTA: {}. Please check that the two files belong to the same genome.",
                missing.join(", ")
            )));
        }
    }

    // 4. run dir + query staging
    let run_dir = state.run_dir(project_id, run_id);
    std::fs::create_dir_all(&run_dir)?;
    let query_ids_json = serde_json::to_string(&query_file_ids).unwrap();
    let _ = query_ids_json;

    ctx.set_step("Preparing the query genomes");
    for (qid, name, path) in &queries {
        let qdir = run_dir.join("queries").join(qid.to_string());
        std::fs::create_dir_all(&qdir)?;
        pipeline::sanitize_into(path, &qdir.join("query.fa"))?;
        ctx.log(&format!("Prepared query genome \u{201c}{name}\u{201d}."));
    }
    if let Some((pname, ppath)) = &panel_path {
        let pdir = run_dir.join("panel");
        std::fs::create_dir_all(&pdir)?;
        pipeline::sanitize_into(ppath, &pdir.join("panel.fa"))?;
        ctx.log(&format!("Prepared gene panel \u{201c}{pname}\u{201d}."));
    }

    // 5. run comparisons in parallel
    let done = Arc::new(AtomicU64::new(0));
    let mut handles = Vec::new();
    let state2 = state.clone();
    let panel_dir = if panel_path.is_some() {
        Some(run_dir.join("panel"))
    } else {
        None
    };
    for (qid, name, _path) in queries.iter() {
        let state3 = state.clone();
        let run_id2 = run_id;
        let project_id2 = project_id;
        let params2 = params.clone();
        let qid = *qid;
        let name = name.clone();
        let name2 = name.clone();
        let ref_fa = staged_fa.clone();
        let ref_gff2 = staged_gff.clone();
        let panel_fa = panel_dir.as_ref().map(|p| p.join("panel.fa"));
        let done2 = done.clone();
        let ctx2 = ctx.clone();
        let handle = tokio::task::spawn_blocking(move || {
            let qdir = state3
                .run_dir(project_id2, run_id2)
                .join("queries")
                .join(qid.to_string());
            let work = qdir.join("work");
            let cache = state3.cache_dir(project_id2);
            let dirs = WorkDirs {
                work: &work,
                cache: &cache,
            };
            let inputs = ComparisonInputs {
                ref_fasta: &ref_fa,
                ref_gff: &ref_gff2,
                qry_fasta: &qdir.join("query.fa"),
                query_name: &name,
                panel_fasta: panel_fa.as_deref(),
                params: &params2,
            };
            let tools = ToolPaths::discover()?;
            let progress = move |msg: &str, _d: u32, _t: u32| {
                ctx2.log(&format!("[{name2}] {msg}"));
                ctx2.set_step(&format!("{msg} ({name2})"));
            };
            let res = pipeline::run_comparison(&tools, &inputs, &dirs, &progress)?;
            Ok((qid, name, res))
        });
        let sem = state2.cpu_slots.clone();
        let handle = async move {
            let permit = sem.acquire_owned().await;
            let r = handle.await;
            drop(permit);
            let d = done2.fetch_add(1, Ordering::SeqCst) + 1;
            let _ = d;
            r
        };
        handles.push(handle);
    }
    let mut results: Vec<(i64, String, ComparisonResult)> = Vec::new();
    for h in handles {
        let outcome: std::result::Result<
            straincompass_engine::Result<(i64, String, ComparisonResult)>,
            tokio::task::JoinError,
        > = h.await;
        match outcome {
            Ok(Ok(triple)) => {
                ctx.log(&format!(
                    "Finished comparing \u{201c}{}\u{201d}: {} genes scored.",
                    triple.2.genes_coverage.len(),
                    triple.1
                ));
                results.push(triple);
            }
            Ok(Err(e)) => return Err(e.into()),
            Err(e) => {
                return Err(crate::error::ApiError::Internal(format!(
                    "A comparison task crashed. ({e})"
                )))
            }
        }
    }
    results.sort_by_key(|(qid, _, _)| *qid);

    // 6. write artifacts
    ctx.set_step("Writing the result tables");
    for (qid, name, res) in &results {
        let qdir = run_dir.join("queries").join(qid.to_string());
        pipeline::write_genes_coverage_tsv(&res.genes_coverage, &qdir.join("genes_coverage.tsv"))?;
        pipeline::write_gaps_tsv(&res.unaligned_gaps, &qdir.join("unaligned_gaps.tsv"))?;
        if let Some(p) = &res.panel {
            pipeline::write_panel_tsv(p, &qdir.join("panel_recheck.tsv"))?;
        }
        if let Some(rep) = &res.dnadiff_report {
            std::fs::write(qdir.join("dnadiff.report"), rep)?;
        }
        std::fs::write(
            qdir.join("result.json"),
            serde_json::to_vec(&ComparisonResultJson {
                query_name: name.clone(),
                genes_coverage: res.genes_coverage.clone(),
                unaligned_gaps: res.unaligned_gaps.clone(),
                panel: res.panel.clone(),
                blocks: res.blocks.clone(),
                ref_lengths: res.ref_lengths.clone(),
            })
            .unwrap(),
        )?;
        // Variant events were precomputed by the comparison itself: cache
        // them now so the alignment viewer never computes on first open.
        std::fs::write(
            qdir.join("variants.json"),
            serde_json::to_vec(&res.events).unwrap(),
        )?;
        ctx.log(&format!(
            "Wrote the result tables for \u{201c}{name}\u{201d} to disk."
        ));
        let _ = name;
    }
    std::fs::write(run_dir.join("params.json"), &params_json)?;
    // reference metadata for viewers
    {
        let genes = straincompass_engine::gff::parse_gff(&staged_gff)?;
        let recs = straincompass_engine::fasta::parse_fasta(&staged_fa)?;
        let lengths: Vec<(String, u64)> = {
            let mut v: Vec<(String, u64)> = recs
                .iter()
                .map(|r| (r.id.clone(), r.seq.len() as u64))
                .collect();
            v.sort();
            v
        };
        let genes_json: Vec<straincompass_types::WgaGene> = genes
            .iter()
            .map(|g| straincompass_types::WgaGene {
                locus_tag: g.locus_tag.clone(),
                symbol: g.symbol.clone(),
                biotype: g.biotype.clone(),
                seqid: g.seqid.clone(),
                start: g.start,
                end: g.end,
                strand: g.strand,
                product: g.product.clone(),
            })
            .collect();
        std::fs::write(
            run_dir.join("reference.json"),
            serde_json::to_vec(&serde_json::json!({
                "lengths": lengths,
                "genes": genes_json,
            }))
            .unwrap(),
        )?;
    }

    // 7. the presence/absence matrix across queries
    if !results.is_empty() {
        let mut rows: Vec<MatrixRow> = Vec::new();
        let n = results.len();
        for (i, g) in results[0].2.genes_coverage.iter().enumerate() {
            let mut calls = Vec::with_capacity(n);
            let mut covs = Vec::with_capacity(n);
            for (_, _, res) in &results {
                let r = &res.genes_coverage[i];
                calls.push(r.call);
                covs.push(r.cov_pct);
            }
            rows.push(MatrixRow {
                locus_tag: g.locus_tag.clone(),
                symbol: g.symbol.clone(),
                biotype: g.biotype.clone(),
                seqid: g.seqid.clone(),
                start: g.start,
                end: g.end,
                calls,
                cov_pcts: covs,
            });
        }
        write_matrix_tsv(
            &rows,
            &results
                .iter()
                .map(|(_, n, _)| n.clone())
                .collect::<Vec<_>>(),
            &run_dir.join("matrix.tsv"),
        )?;
        ctx.log("Built the presence/absence table across all queries.");
    }

    // 8. done
    {
        let conn = state.db.lock().unwrap();
        conn.execute(
            "UPDATE runs SET status = 'succeeded', step = NULL,
             finished_at = strftime('%Y-%m-%dT%H:%M:%fZ','now') WHERE id = ?1",
            [run_id],
        )?;
    }
    ctx.log("Run finished successfully.");
    Ok(())
}

#[derive(serde::Serialize, serde::Deserialize)]
pub struct ComparisonResultJson {
    pub query_name: String,
    pub genes_coverage: Vec<straincompass_types::GeneCoverageRow>,
    pub unaligned_gaps: Vec<straincompass_types::GapRow>,
    pub panel: Option<Vec<straincompass_types::PanelRow>>,
    pub blocks: Vec<straincompass_types::WgaBlock>,
    pub ref_lengths: Vec<(String, u64)>,
}

fn write_matrix_tsv(
    rows: &[MatrixRow],
    query_names: &[String],
    out: &std::path::Path,
) -> std::io::Result<()> {
    let mut w = std::io::BufWriter::new(std::fs::File::create(out)?);
    write!(w, "locus_tag\tsymbol\tbiotype\tseqid\tstart\tend")?;
    for n in query_names {
        write!(w, "\t{n}")?;
    }
    writeln!(w)?;
    for r in rows {
        write!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}",
            r.locus_tag, r.symbol, r.biotype, r.seqid, r.start, r.end
        )?;
        for c in &r.calls {
            write!(w, "\t{}", c.as_str())?;
        }
        writeln!(w)?;
    }
    Ok(())
}

/// Load the per-query result of a finished run (for the result endpoints).
pub fn load_query_result(
    state: &SharedState,
    project_id: i64,
    run_id: i64,
    query_file_id: i64,
) -> ApiResult<ComparisonResultJson> {
    let path = state
        .run_dir(project_id, run_id)
        .join("queries")
        .join(query_file_id.to_string())
        .join("result.json");
    let bytes = std::fs::read(&path).map_err(|_| {
        crate::error::ApiError::NotFound(
            "The results for this query are not available (the run may not have finished).".into(),
        )
    })?;
    let mut res: ComparisonResultJson = serde_json::from_slice(&bytes)
        .map_err(|_| crate::error::ApiError::Internal("A result file is unreadable.".into()))?;
    backfill_protein_ids(state, project_id, &mut res.genes_coverage);
    Ok(res)
}

/// Slim read of a finished query: just the display name and the
/// alignment blocks. The alignment viewer does not need the big
/// gene-coverage rows, so this skips deserializing them (and the
/// protein-id GFF backfill `load_query_result` sometimes triggers).
pub fn load_query_meta(
    state: &SharedState,
    project_id: i64,
    run_id: i64,
    query_file_id: i64,
) -> ApiResult<(String, Vec<straincompass_types::WgaBlock>)> {
    let path = state
        .run_dir(project_id, run_id)
        .join("queries")
        .join(query_file_id.to_string())
        .join("result.json");
    let bytes = std::fs::read(&path).map_err(|_| {
        crate::error::ApiError::NotFound(
            "The results for this query are not available (the run may not have finished).".into(),
        )
    })?;
    #[derive(serde::Deserialize)]
    struct Meta {
        query_name: String,
        #[serde(default)]
        blocks: Vec<straincompass_types::WgaBlock>,
    }
    let meta: Meta = serde_json::from_slice(&bytes)
        .map_err(|_| crate::error::ApiError::Internal("A result file is unreadable.".into()))?;
    Ok((meta.query_name, meta.blocks))
}

/// Runs computed before protein accessions were persisted carry empty
/// `protein_id`. Join them in from the project's reference annotation at read
/// time, keyed by locus tag (then `old_locus_tag`, for re-annotated GFFs), so
/// those runs show the Protein column too. New runs already persist it.
fn backfill_protein_ids(
    state: &SharedState,
    project_id: i64,
    rows: &mut [straincompass_types::GeneCoverageRow],
) {
    if rows.is_empty() || rows.iter().all(|r| !r.protein_id.is_empty()) {
        return;
    }
    let staged = state
        .project_dir(project_id)
        .join("reference")
        .join("ref.gff");
    let gff_path = match std::fs::exists(&staged) {
        Ok(true) => staged,
        _ => match crate::routes::uploads::reference_paths(state, project_id) {
            Ok(Some((_, gff))) => gff,
            _ => {
                tracing::warn!("no reference GFF to backfill protein accessions");
                return;
            }
        },
    };
    let genes = match straincompass_engine::gff::parse_gff(&gff_path) {
        Ok(g) => g,
        Err(e) => {
            tracing::warn!("cannot parse reference GFF for protein backfill: {e}");
            return;
        }
    };
    let by_tag: std::collections::HashMap<&str, &str> = genes
        .iter()
        .filter(|g| !g.protein_id.is_empty())
        .map(|g| (g.locus_tag.as_str(), g.protein_id.as_str()))
        .collect();
    let by_old_tag: std::collections::HashMap<&str, &str> = genes
        .iter()
        .filter(|g| !g.protein_id.is_empty() && !g.old_locus_tag.is_empty())
        .map(|g| (g.old_locus_tag.as_str(), g.protein_id.as_str()))
        .collect();
    for r in rows.iter_mut() {
        if r.protein_id.is_empty() {
            r.protein_id = by_tag
                .get(r.locus_tag.as_str())
                .or_else(|| by_old_tag.get(r.locus_tag.as_str()))
                .map(|p| (*p).to_string())
                .unwrap_or_default();
        }
    }
}

pub fn load_reference_json(
    state: &SharedState,
    project_id: i64,
    run_id: i64,
) -> ApiResult<serde_json::Value> {
    let path = state.run_dir(project_id, run_id).join("reference.json");
    let bytes = std::fs::read(&path).map_err(|_| {
        crate::error::ApiError::NotFound("The reference data for this run is not available.".into())
    })?;
    serde_json::from_slice(&bytes)
        .map_err(|_| crate::error::ApiError::Internal("A result file is unreadable.".into()))
}

/// Resolve (project_id, query_file_ids, status) for a run.
pub fn run_meta(state: &SharedState, run_id: i64) -> ApiResult<(i64, Vec<i64>, String)> {
    let conn = state.db.lock().unwrap();
    conn.query_row(
        "SELECT project_id, query_ids, status FROM runs WHERE id = ?1",
        [run_id],
        |r| {
            Ok((
                r.get::<_, i64>(0)?,
                r.get::<_, String>(1)?,
                r.get::<_, String>(2)?,
            ))
        },
    )
    .map(|(p, q, s)| (p, serde_json::from_str(&q).unwrap_or_default(), s))
    .map_err(|_| crate::error::ApiError::NotFound("This run does not exist (anymore).".into()))
}

/// Build MSA sources for gene_detail.
pub fn msa_sources(
    state: &SharedState,
    project_id: i64,
    run_id: i64,
) -> ApiResult<(PathBuf, PathBuf, RunParams, Vec<QueryAlignmentSourceOwned>)> {
    let (_, query_ids, _) = run_meta(state, run_id)?;
    let run_dir = state.run_dir(project_id, run_id);
    let params: RunParams = {
        let conn = state.db.lock().unwrap();
        let json: String = conn.query_row(
            "SELECT params_json FROM runs WHERE id = ?1",
            [run_id],
            |r| r.get(0),
        )?;
        serde_json::from_str(&json).map_err(|_| {
            crate::error::ApiError::Internal("The stored parameters are unreadable.".into())
        })?
    };
    let mut sources = Vec::new();
    for qid in query_ids {
        let name: Option<String> = {
            let conn = state.db.lock().unwrap();
            conn.query_row("SELECT display_name FROM files WHERE id = ?1", [qid], |r| {
                r.get(0)
            })
            .ok()
        };
        let Some(name) = name else { continue };
        let qdir = run_dir.join("queries").join(qid.to_string());
        let qry = qdir.join("query.fa");
        let delta = qdir.join("work").join("cmp.delta");
        if qry.is_file() && delta.is_file() {
            sources.push(QueryAlignmentSourceOwned {
                query_id: qid,
                query_name: name,
                qry_fasta: qry,
                delta,
            });
        }
    }
    let ref_dir = state.project_dir(project_id).join("reference");
    Ok((
        ref_dir.join("ref.fa"),
        ref_dir.join("ref.gff"),
        params,
        sources,
    ))
}

pub struct QueryAlignmentSourceOwned {
    pub query_id: i64,
    pub query_name: String,
    pub qry_fasta: PathBuf,
    pub delta: PathBuf,
}

impl QueryAlignmentSourceOwned {
    pub fn borrow(&self) -> QueryAlignmentSource<'_> {
        QueryAlignmentSource {
            query_id: self.query_id,
            query_name: self.query_name.clone(),
            qry_fasta: &self.qry_fasta,
            delta: &self.delta,
        }
    }
}

/// Per-query guards for the on-demand variant event computation, keyed
/// by the cache file path: two requests for the SAME query do not
/// duplicate the heavy walk, while different queries compute in
/// parallel (the alignment viewer fans out one blocking task per query).
static VARIANT_LOCKS: std::sync::OnceLock<
    std::sync::Mutex<std::collections::HashMap<PathBuf, std::sync::Arc<std::sync::Mutex<()>>>>,
> = std::sync::OnceLock::new();

fn variant_lock(path: &std::path::Path) -> std::sync::Arc<std::sync::Mutex<()>> {
    let map = VARIANT_LOCKS.get_or_init(|| std::sync::Mutex::new(std::collections::HashMap::new()));
    let mut guard = map.lock().unwrap();
    std::sync::Arc::clone(guard.entry(path.to_path_buf()).or_default())
}

/// Load the per-base variant events of one query against the reference.
/// Computed on demand from the kept delta file on first use, then cached
/// as variants.json next to the query's result.json. Queries whose
/// working files are gone get empty events (the viewer still shows
/// their alignment blocks). Runs now write variants.json on completion,
/// so this on-demand path is a backfill for runs made before that.
///
/// `ref_records` comes pre-parsed from the caller, which batches one of
/// these calls per query in parallel: the multi-megabase reference
/// fasta is parsed once for the whole run instead of once per query.
pub fn load_variants(
    state: &SharedState,
    project_id: i64,
    run_id: i64,
    query_file_id: i64,
    ref_records: &[straincompass_engine::fasta::FastaRecord],
) -> ApiResult<std::collections::BTreeMap<String, straincompass_types::AlignmentEvents>> {
    let qdir = state
        .run_dir(project_id, run_id)
        .join("queries")
        .join(query_file_id.to_string());
    let path = qdir.join("variants.json");
    if !path.is_file() {
        let guard = variant_lock(&path);
        let _guard = guard.lock().unwrap();
        if !path.is_file() {
            let qry_fa = qdir.join("query.fa");
            let delta = qdir.join("work").join("cmp.delta");
            if qry_fa.is_file() && delta.is_file() {
                straincompass_engine::variants::write_variant_events_with_ref(
                    ref_records,
                    &qry_fa,
                    &delta,
                    &path,
                )?;
            } else {
                return Ok(Default::default());
            }
        }
    }
    let bytes = std::fs::read(&path).map_err(|_| {
        crate::error::ApiError::Internal(
            "The variant data for this query could not be read.".into(),
        )
    })?;
    serde_json::from_slice(&bytes).map_err(|_| {
        crate::error::ApiError::Internal("The variant data for this query is unreadable.".into())
    })
}
