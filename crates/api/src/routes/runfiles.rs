//! Run artifacts surfaced through the API ("Files" panel) and the gene
//! alignment export (FASTA / clustal).

use crate::error::{ApiError, ApiResult};
use crate::jobs;
use crate::state::SharedState;
use axum::body::Body;
use axum::extract::{Path, Query, State};
use axum::http::header;
use axum::http::HeaderMap;
use axum::response::{IntoResponse, Response};
use serde::Deserialize;
use serde_json::json;
use std::path::PathBuf;

/// GET /runs/{id}/files : artifact list with friendly names.
pub async fn list_run_files(
    State(state): State<SharedState>,
    Path(run_id): Path<i64>,
) -> ApiResult<axum::Json<serde_json::Value>> {
    let (project_id, query_ids, status) = jobs::run_meta(&state, run_id)?;
    let run_dir = state.run_dir(project_id, run_id);
    let mut files = Vec::new();
    // parameters + log
    push_artifact(&mut files, &run_dir.join("params.json"), "Parameters used (JSON)", "parameters.json");
    // per query artifacts
    let conn = state.db.lock().unwrap();
    let mut names = Vec::new();
    for qid in &query_ids {
        let name: Option<String> = conn
            .query_row("SELECT display_name FROM files WHERE id = ?1", [qid], |r| r.get(0))
            .ok();
        names.push((qid, name));
    }
    drop(conn);
    for (qid, name) in names {
        let Some(qname) = name else { continue };
        let qdir = run_dir.join("queries").join(qid.to_string());
        let stem = sanitize_stem(&qname);
        push_artifact(&mut files, &qdir.join("genes_coverage.tsv"), "Genes coverage table", &format!("{stem}_genes_coverage.tsv"));
        push_artifact(&mut files, &qdir.join("unaligned_gaps.tsv"), "Unaligned gaps table", &format!("{stem}_unaligned_gaps.tsv"));
        push_artifact(&mut files, &qdir.join("panel_recheck.tsv"), "Gene panel recheck table", &format!("{stem}_panel_recheck.tsv"));
        push_artifact(&mut files, &qdir.join("dnadiff.report"), "Overall alignment report", &format!("{stem}_dnadiff.report"));
    }
    if query_ids.len() >= 1 {
        push_artifact(&mut files, &run_dir.join("matrix.tsv"), "Presence/absence table (all queries)", "presence_absence_matrix.tsv");
    }
    Ok(axum::Json(json!({
        "run_id": run_id,
        "status": status,
        "files": files,
    })))
}

fn sanitize_stem(name: &str) -> String {
    let stem = name.trim_end_matches(".fasta").trim_end_matches(".fa").trim_end_matches(".fna");
    stem.chars()
        .map(|c| if c.is_ascii_alphanumeric() || matches!(c, '.' | '-' | '_') { c } else { '_' })
        .collect()
}

fn push_artifact(out: &mut Vec<serde_json::Value>, path: &std::path::Path, friendly: &str, name: &str) {
    if let Ok(meta) = std::fs::metadata(path) {
        if meta.is_file() {
            out.push(json!({
                "name": name,
                "friendly": friendly,
                "size": meta.len(),
                "human_size": crate::files::human_size(meta.len()),
            }));
        }
    }
}

/// GET /runs/{id}/files/{name} : view or download one artifact.
pub async fn get_run_file(
    State(state): State<SharedState>,
    Path((run_id, name)): Path<(i64, String)>,
    Query(dl): Query<DownloadQ>,
) -> ApiResult<Response> {
    let (project_id, query_ids, _) = jobs::run_meta(&state, run_id)?;
    let run_dir = state.run_dir(project_id, run_id);
    let candidates = candidate_paths(&run_dir, &query_ids);
    let Some((path, fname)) = candidates.iter().find(|(_, n)| n == &name) else {
        return Err(ApiError::NotFound("This file does not exist.".into()));
    };
    let path: PathBuf = path.clone();
    let bytes = std::fs::read(&path)
        .map_err(|_| ApiError::NotFound("This file does not exist.".into()))?;
    let mut headers = HeaderMap::new();
    let mime = if fname.ends_with(".json") {
        "application/json"
    } else if fname.ends_with(".tsv") {
        "text/tab-separated-values; charset=utf-8"
    } else if fname.ends_with(".csv") {
        "text/csv; charset=utf-8"
    } else {
        "text/plain; charset=utf-8"
    };
    headers.insert(header::CONTENT_TYPE, mime.parse().unwrap());
    if dl.download.unwrap_or(false) {
        headers.insert(
            header::CONTENT_DISPOSITION,
            format!("attachment; filename=\"{fname}\"").parse().unwrap(),
        );
    }
    Ok((headers, Body::from(bytes)).into_response())
}

#[derive(Deserialize)]
pub struct DownloadQ {
    pub download: Option<bool>,
}

fn candidate_paths(run_dir: &std::path::Path, query_ids: &[i64]) -> Vec<(PathBuf, String)> {
    let mut out = vec![(run_dir.join("params.json"), "parameters.json".to_string())];
    // Query names are unknown here; mirror the list_run_files naming by
    // scanning the query dirs.
    if let Ok(conn_placeholder) = scan_dirs(run_dir, query_ids) {
        out.extend(conn_placeholder);
    }
    out.push((
        run_dir.join("matrix.tsv"),
        "presence_absence_matrix.tsv".to_string(),
    ));
    out
}

fn scan_dirs(run_dir: &std::path::Path, query_ids: &[i64]) -> ApiResult<Vec<(PathBuf, String)>> {
    let mut out = Vec::new();
    let queries_dir = run_dir.join("queries");
    if let Ok(entries) = std::fs::read_dir(&queries_dir) {
        for e in entries.flatten() {
            let qdir = e.path();
            let qid = e
                .file_name()
                .to_string_lossy()
                .to_string();
            if !query_ids.iter().any(|q| q.to_string() == qid) {
                continue;
            }
            // read the query name from result.json
            let name = std::fs::read(qdir.join("result.json"))
                .ok()
                .and_then(|b| serde_json::from_slice::<serde_json::Value>(&b).ok())
                .and_then(|v| v["query_name"].as_str().map(|s| s.to_string()))
                .unwrap_or_else(|| format!("query_{qid}"));
            let stem = sanitize_stem(&name);
            out.push((qdir.join("genes_coverage.tsv"), format!("{stem}_genes_coverage.tsv")));
            out.push((qdir.join("unaligned_gaps.tsv"), format!("{stem}_unaligned_gaps.tsv")));
            out.push((qdir.join("panel_recheck.tsv"), format!("{stem}_panel_recheck.tsv")));
            out.push((qdir.join("dnadiff.report"), format!("{stem}_dnadiff.report")));
        }
    }
    Ok(out)
}

/// GET /runs/{id}/gene/{locus}/export?format=fasta|clustal
pub async fn export_gene_alignment(
    State(state): State<SharedState>,
    Path((run_id, locus)): Path<(i64, String)>,
    Query(q): Query<GeneExportQ>,
) -> ApiResult<Response> {
    let (project_id, _, _) = jobs::run_meta(&state, run_id)?;
    let (ref_fa, ref_gff, params, sources) = jobs::msa_sources(&state, project_id, run_id)?;
    let borrowed: Vec<_> = sources.iter().map(|s| s.borrow()).collect();
    let detail = bactiment_engine::pipeline::gene_detail(&ref_fa, &ref_gff, &params, &locus, &borrowed)?;
    let clustal = q.format.as_deref() == Some("clustal");
    let (body, ext, mime) = if clustal {
        (clustal_format(&detail), "aln", "text/plain; charset=utf-8")
    } else {
        (fasta_format(&detail), "fasta", "text/plain; charset=utf-8")
    };
    let mut headers = HeaderMap::new();
    headers.insert(header::CONTENT_TYPE, mime.parse().unwrap());
    headers.insert(
        header::CONTENT_DISPOSITION,
        format!(
            "attachment; filename=\"{}_gene_alignment.{}\"",
            sanitize_stem(&locus),
            ext
        )
        .parse()
        .unwrap(),
    );
    Ok((headers, Body::from(body)).into_response())
}

#[derive(Deserialize)]
pub struct GeneExportQ {
    pub format: Option<String>,
}

fn fasta_format(d: &bactiment_types::GeneDetail) -> String {
    let mut out = format!(">{}_reference\n", d.locus_tag);
    for q in &d.queries {
        out.push_str(&format!(
            ">{}__{} {}\n",
            d.locus_tag, q.query_name, q.query_name
        ));
    }
    // Build per-query concatenated alignment (reference row + query rows).
    // Walk block by block; queries missing a block get gaps.
    let mut ref_row: Vec<u8> = Vec::new();
    let mut qry_rows: Vec<Vec<u8>> = vec![Vec::new(); d.queries.len()];
    let mut offsets = vec![0usize; d.queries.len()];
    for qi in 0..d.queries.len() {
        offsets[qi] = 0;
    }
    // Collect all block boundaries in reference order across queries:
    // simpler approach: concatenate each query's blocks independently,
    // since blocks are per query.
    for (qi, q) in d.queries.iter().enumerate() {
        if qi == 0 {
            for b in &q.blocks {
                ref_row.extend_from_slice(b.ref_seq.as_bytes());
            }
        }
        // reference pieces may differ between queries when coverage
        // differs; pad is impossible without global coordinates, so use
        // each query's own reference slice.
        for b in &q.blocks {
            qry_rows[qi].extend_from_slice(b.qry_seq.as_bytes());
        }
    }
    // reference sequence from the first query that has blocks
    let mut out = String::new();
    let mut label = format!("{}_reference", d.locus_tag);
    out.push_str(&format!(">{label}\n"));
    for chunk in ref_row.chunks(60) {
        out.push_str(&String::from_utf8_lossy(chunk));
        out.push('\n');
    }
    for (qi, q) in d.queries.iter().enumerate() {
        if qry_rows[qi].is_empty() {
            continue;
        }
        label = format!("{}__{}", d.locus_tag, q.query_name);
        out.push_str(&format!(">{label}\n"));
        for chunk in qry_rows[qi].chunks(60) {
            out.push_str(&String::from_utf8_lossy(chunk));
            out.push('\n');
        }
    }
    out
}

fn clustal_format(d: &bactiment_types::GeneDetail) -> String {
    // simple clustal-like output: one block per alignment block, aligned
    // rows per query, 60 columns per line
    let mut out = String::from("CLUSTAL W (bactiment) multiple sequence alignment\n\n");
    let mut block_no = 0;
    let max_name = d
        .queries
        .iter()
        .map(|q| q.query_name.len())
        .chain(std::iter::once(8))
        .max()
        .unwrap_or(8)
        + 4;
    // for each query, print reference vs query rows per block
    for q in &d.queries {
        for b in &q.blocks {
            block_no += 1;
            out.push_str(&format!("Alignment block {} (reference {}-{}, identity {:.1}%)\n", block_no, b.ref_start, b.ref_end, b.identity));
            let rows = [
                (format!("{}_reference", d.locus_tag), b.ref_seq.clone()),
                (q.query_name.clone(), b.qry_seq.clone()),
            ];
            let mut pos = 0usize;
            while pos < b.ref_seq.len() {
                let end = (pos + 60).min(b.ref_seq.len());
                for (name, seq) in &rows {
                    out.push_str(&format!("{:<width$}{}\n", name, &seq[pos..end], width = max_name));
                }
                out.push('\n');
                pos = end;
            }
        }
    }
    out
}
