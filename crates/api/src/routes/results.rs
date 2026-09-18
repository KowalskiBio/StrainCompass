//! Result table endpoints with server side sort/filter/pagination,
//! the WGA viewer data and the gene alignment (MSA) endpoint.

use crate::error::{ApiError, ApiResult};
use crate::jobs;
use crate::state::SharedState;
use axum::extract::{Path, Query, State};
use axum::Json;
use serde_json::Value;
use straincompass_types::{
    Call, GainedAnchor, GainedRow, GapRow, GeneCoverageRow, GeneDetail, MatrixRow, Page, PanelRow,
    TableQuery, WgaData, WgaQuery,
};

fn resolve_query_id(
    state: &SharedState,
    run_id: i64,
    query_id: Option<i64>,
) -> ApiResult<(i64, i64, i64)> {
    let (project_id, query_ids, status) = jobs::run_meta(state, run_id)?;
    if status != "succeeded" {
        return Err(ApiError::BadRequest(
            "This run has not finished yet. Please wait for it to complete.".into(),
        ));
    }
    let qid = match query_id {
        Some(q) => {
            if !query_ids.contains(&q) {
                return Err(ApiError::BadRequest(
                    "This query is not part of the run.".into(),
                ));
            }
            q
        }
        None => *query_ids
            .first()
            .ok_or_else(|| ApiError::BadRequest("This run has no queries.".into()))?,
    };
    Ok((project_id, run_id, qid))
}

fn matches_search(row_text: &str, search: &Option<String>) -> bool {
    match search {
        Some(s) if !s.trim().is_empty() => {
            row_text.to_lowercase().contains(&s.trim().to_lowercase())
        }
        _ => true,
    }
}

fn call_filter_ok(call: Call, filter: &Option<String>) -> bool {
    match filter.as_deref() {
        Some("present") => call == Call::Present,
        Some("partial") => call == Call::Partial,
        Some("absent") => call == Call::Absent,
        _ => true,
    }
}

fn row_str(row: &GeneCoverageRow) -> String {
    format!("{} {} {}", row.locus_tag, row.symbol, row.seqid)
}

fn sort_rows(rows: &mut [GeneCoverageRow], by: &Option<String>, dir: &Option<String>) {
    let asc = dir.as_deref() != Some("desc");
    let by = by.clone().unwrap_or_default();
    rows.sort_by(|a, b| {
        let ord = match by.as_str() {
            "locus_tag" => a.locus_tag.cmp(&b.locus_tag),
            "symbol" => a.symbol.cmp(&b.symbol),
            "biotype" => a.biotype.cmp(&b.biotype),
            "seqid" => a.seqid.cmp(&b.seqid),
            "start" => a.start.cmp(&b.start),
            "end" => a.end.cmp(&b.end),
            "length" => a.length.cmp(&b.length),
            "cov_bp" => a.cov_bp.cmp(&b.cov_bp),
            "cov_pct" => a
                .cov_pct
                .partial_cmp(&b.cov_pct)
                .unwrap_or(std::cmp::Ordering::Equal),
            "call" => a.call.as_str().cmp(b.call.as_str()),
            "best_identity" => a
                .best_identity
                .partial_cmp(&b.best_identity)
                .unwrap_or(std::cmp::Ordering::Equal),
            "mismatches" => a.mismatches.cmp(&b.mismatches),
            "indels" => a.indels.cmp(&b.indels),
            _ => std::cmp::Ordering::Equal,
        };
        if asc {
            ord
        } else {
            ord.reverse()
        }
    });
}

/// GET /runs/{id}/genes_coverage
pub async fn genes_coverage(
    State(state): State<SharedState>,
    Path(run_id): Path<i64>,
    Query(q): Query<TableQuery>,
) -> ApiResult<Json<Page<GeneCoverageRow>>> {
    let (project_id, run_id, qid) = resolve_query_id(&state, run_id, q.query_id)?;
    let res = jobs::load_query_result(&state, project_id, run_id, qid)?;
    let mut rows = res.genes_coverage;
    rows.retain(|r| matches_search(&row_str(r), &q.search) && call_filter_ok(r.call, &q.call));
    sort_rows(&mut rows, &q.sort_by, &q.sort_dir);
    let total = rows.len() as u64;
    let page = apply_page(&rows, q.page, q.page_size);
    Ok(Json(Page { rows: page, total }))
}

fn apply_page<T: Clone>(rows: &[T], page: u64, page_size: u64) -> Vec<T> {
    let page_size = page_size.clamp(1, 1000);
    let start = (page * page_size) as usize;
    rows.iter()
        .skip(start)
        .take(page_size as usize)
        .cloned()
        .collect()
}

/// GET /runs/{id}/unaligned_gaps
pub async fn unaligned_gaps(
    State(state): State<SharedState>,
    Path(run_id): Path<i64>,
    Query(q): Query<TableQuery>,
) -> ApiResult<Json<Page<GapRow>>> {
    let (project_id, run_id, qid) = resolve_query_id(&state, run_id, q.query_id)?;
    let res = jobs::load_query_result(&state, project_id, run_id, qid)?;
    let mut rows = res.unaligned_gaps;
    rows.retain(|r| matches_search(&format!("{} {}", r.seqid, r.genes.join(" ")), &q.search));
    let asc = q.sort_dir.as_deref() != Some("desc");
    let by = q.sort_by.clone().unwrap_or_default();
    rows.sort_by(|a, b| {
        let ord = match by.as_str() {
            "seqid" => a.seqid.cmp(&b.seqid),
            "length" => a.length.cmp(&b.length),
            "start" => a.start.cmp(&b.start),
            "end" => a.end.cmp(&b.end),
            "n_genes" => a.n_genes.cmp(&b.n_genes),
            _ => (a.seqid.clone(), a.start).cmp(&(b.seqid.clone(), b.start)),
        };
        if asc {
            ord
        } else {
            ord.reverse()
        }
    });
    let total = rows.len() as u64;
    let page = apply_page(&rows, q.page, q.page_size);
    Ok(Json(Page { rows: page, total }))
}

/// GET /runs/{id}/gained : stretches of one query with no alignment to the
/// reference, and the genes predicted inside them.
pub async fn gained(
    State(state): State<SharedState>,
    Path(run_id): Path<i64>,
    Query(q): Query<TableQuery>,
) -> ApiResult<Json<Page<GainedRow>>> {
    let (project_id, run_id, qid) = resolve_query_id(&state, run_id, q.query_id)?;
    let res = jobs::load_query_result(&state, project_id, run_id, qid)?;
    // A run from before this feature has no gained rows and cannot grow
    // them without being re-run: the gene finder may be missing and the old
    // params never recorded a minimum length, so a backfill would cache a
    // result the user never asked for.
    let mut rows = res.gained.ok_or_else(|| {
        ApiError::BadRequest(
            "This run was computed before gained regions were available. Please run the comparison again to see them."
                .into(),
        )
    })?;
    rows.retain(|r| {
        matches_search(
            &format!(
                "{} {} {} {}",
                r.qry_seqid, r.anchor_seqid, r.left_gene, r.right_gene
            ),
            &q.search,
        )
    });
    // The anchor filter reuses `call`, as matrix() already does for
    // "not_present", rather than adding a field to the query struct the
    // export endpoints share.
    match q.call.as_deref() {
        Some("anchored") => rows.retain(|r| r.anchor != GainedAnchor::Unanchored),
        Some("unanchored") => rows.retain(|r| r.anchor == GainedAnchor::Unanchored),
        _ => {}
    }
    let asc = q.sort_dir.as_deref() != Some("desc");
    let by = q.sort_by.clone().unwrap_or_default();
    rows.sort_by(|a, b| {
        let ord = match by.as_str() {
            "start" => a.start.cmp(&b.start),
            "end" => a.end.cmp(&b.end),
            "length" => a.length.cmp(&b.length),
            "gc_pct" => a.gc_pct.total_cmp(&b.gc_pct),
            "anchor" => a.anchor.as_str().cmp(b.anchor.as_str()),
            "anchor_seqid" => a.anchor_seqid.cmp(&b.anchor_seqid),
            "anchor_start" => a.anchor_start.cmp(&b.anchor_start),
            // An absent count sorts last whichever way the column is
            // pointing: "not measured" is not a small number.
            "n_orfs" => none_last(a.n_orfs, b.n_orfs, asc),
            "n_orfs_complete" => none_last(a.n_orfs_complete, b.n_orfs_complete, asc),
            _ => (a.qry_seqid.clone(), a.start).cmp(&(b.qry_seqid.clone(), b.start)),
        };
        if asc {
            ord
        } else {
            ord.reverse()
        }
    });
    let total = rows.len() as u64;
    let page = apply_page(&rows, q.page, q.page_size);
    Ok(Json(Page { rows: page, total }))
}

/// Order two optional counts so that `None` ends up at the end of the table
/// regardless of direction. The caller reverses the result for a descending
/// sort, so the comparison is pre-flipped here.
fn none_last(a: Option<u32>, b: Option<u32>, asc: bool) -> std::cmp::Ordering {
    use std::cmp::Ordering;
    match (a, b) {
        (Some(x), Some(y)) => x.cmp(&y),
        (Some(_), None) => {
            if asc {
                Ordering::Less
            } else {
                Ordering::Greater
            }
        }
        (None, Some(_)) => {
            if asc {
                Ordering::Greater
            } else {
                Ordering::Less
            }
        }
        (None, None) => Ordering::Equal,
    }
}

#[derive(serde::Deserialize)]
pub struct GainedVerifyQuery {
    pub query_id: Option<i64>,
    pub seqid: String,
    pub start: u64,
    pub end: u64,
}

/// GET /runs/{id}/gained/verify?query_id&seqid&start&end : blastn one
/// gained region back against the reference genome.
///
/// "No alignment to the reference" is the anchor-based definition of a
/// gained region and is weaker than "absent from the reference"; this is
/// the closer test, run on demand because it costs a blast database and
/// a search per region and most rows are never questioned.
pub async fn gained_verify(
    State(state): State<SharedState>,
    Path(run_id): Path<i64>,
    Query(q): Query<GainedVerifyQuery>,
) -> ApiResult<Json<straincompass_types::GainedVerify>> {
    let (project_id, run_id, qid) = resolve_query_id(&state, run_id, q.query_id)?;
    let res = jobs::load_query_result(&state, project_id, run_id, qid)?;
    let rows = res.gained.ok_or_else(|| {
        ApiError::BadRequest(
            "This run was computed before gained regions were available. Please run the comparison again to see them."
                .into(),
        )
    })?;
    // The endpoint reads the query fasta, so it answers only for regions
    // this run actually reported - not for arbitrary coordinates.
    if !rows
        .iter()
        .any(|r| r.qry_seqid == q.seqid && r.start == q.start && r.end == q.end)
    {
        return Err(ApiError::BadRequest(
            "This region is not in this query's gained table.".into(),
        ));
    }

    let run_dir = state.run_dir(project_id, run_id);
    let qry_fa = run_dir
        .join("queries")
        .join(qid.to_string())
        .join("query.fa");
    let ref_fa = state
        .project_dir(project_id)
        .join("reference")
        .join("ref.fa");
    for f in [&qry_fa, &ref_fa] {
        if !f.is_file() {
            return Err(ApiError::NotFound(
                "This run's input files are no longer on the server, so the region cannot be re-examined.".into(),
            ));
        }
    }

    let work = run_dir
        .join("queries")
        .join(qid.to_string())
        .join("work")
        .join(format!("gained_verify_{}", uuid::Uuid::new_v4()));
    let scratch = work.clone();
    let seqid = q.seqid;
    let (start, end) = (q.start, q.end);
    // One check at a time holds a cpu slot, like the alignment backfill:
    // the blast search is seconds, but it is still a whole-genome search.
    let _permit = state
        .cpu_slots
        .clone()
        .acquire_owned()
        .await
        .map_err(|_| ApiError::Internal("The server is shutting down.".into()))?;
    let result = tokio::task::spawn_blocking(move || {
        let tools = straincompass_engine::tools::ToolPaths::discover()?;
        straincompass_engine::blast::gained_verify(
            &tools, &ref_fa, &qry_fa, &seqid, start, end, &work,
        )
        .map_err(ApiError::from)
    })
    .await
    .map_err(|e| ApiError::Internal(format!("The reference check crashed. ({e})")))?;
    drop(_permit);
    // The scratch dir is per call; leaving it behind would grow the run
    // by one blast database per click.
    let _ = std::fs::remove_dir_all(&scratch);
    result.map(Json)
}

/// GET /runs/{id}/panel_recheck
pub async fn panel_recheck(
    State(state): State<SharedState>,
    Path(run_id): Path<i64>,
    Query(q): Query<TableQuery>,
) -> ApiResult<Json<Page<PanelRow>>> {
    let (project_id, run_id, qid) = resolve_query_id(&state, run_id, q.query_id)?;
    let res = jobs::load_query_result(&state, project_id, run_id, qid)?;
    let mut rows = res.panel.ok_or_else(|| {
        ApiError::BadRequest("This run has no gene panel results (no panel was provided).".into())
    })?;
    rows.retain(|r| matches_search(&r.gene_id, &q.search) && call_filter_ok(r.call, &q.call));
    let asc = q.sort_dir.as_deref() != Some("desc");
    let by = q.sort_by.clone().unwrap_or_default();
    rows.sort_by(|a, b| {
        let ord = match by.as_str() {
            "gene_id" => a.gene_id.cmp(&b.gene_id),
            "qlen" => a.qlen.cmp(&b.qlen),
            "cov_pct" => a
                .cov_pct
                .partial_cmp(&b.cov_pct)
                .unwrap_or(std::cmp::Ordering::Equal),
            "identity" => a
                .identity
                .partial_cmp(&b.identity)
                .unwrap_or(std::cmp::Ordering::Equal),
            "call" => a.call.as_str().cmp(b.call.as_str()),
            _ => a.gene_id.cmp(&b.gene_id),
        };
        if asc {
            ord
        } else {
            ord.reverse()
        }
    });
    let total = rows.len() as u64;
    let page = apply_page(&rows, q.page, q.page_size);
    Ok(Json(Page { rows: page, total }))
}

/// GET /runs/{id}/matrix : presence/absence across all queries.
pub async fn matrix(
    State(state): State<SharedState>,
    Path(run_id): Path<i64>,
    Query(q): Query<TableQuery>,
) -> ApiResult<Json<Page<MatrixRow>>> {
    let (project_id, query_ids, status) = jobs::run_meta(&state, run_id)?;
    if status != "succeeded" {
        return Err(ApiError::BadRequest(
            "This run has not finished yet. Please wait for it to complete.".into(),
        ));
    }
    let mut results = Vec::new();
    for qid in &query_ids {
        results.push(jobs::load_query_result(&state, project_id, run_id, *qid)?);
    }
    let n = results.len();
    let mut rows: Vec<MatrixRow> = Vec::new();
    if let Some(first) = results.first() {
        for (i, g) in first.genes_coverage.iter().enumerate() {
            let mut calls = Vec::with_capacity(n);
            let mut covs = Vec::with_capacity(n);
            for res in &results {
                calls.push(res.genes_coverage[i].call);
                covs.push(res.genes_coverage[i].cov_pct);
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
    }
    rows.retain(|r| {
        matches_search(
            &format!("{} {} {}", r.locus_tag, r.symbol, r.seqid),
            &q.search,
        )
    });
    if q.call.as_deref() == Some("not_present") {
        rows.retain(|r| r.calls.iter().any(|c| *c != Call::Present));
    }
    let asc = q.sort_dir.as_deref() != Some("desc");
    let by = q.sort_by.clone().unwrap_or_default();
    rows.sort_by(|a, b| {
        let ord = match by.as_str() {
            "symbol" => a.symbol.cmp(&b.symbol),
            "biotype" => a.biotype.cmp(&b.biotype),
            _ => a.locus_tag.cmp(&b.locus_tag),
        };
        if asc {
            ord
        } else {
            ord.reverse()
        }
    });
    let total = rows.len() as u64;
    let page = apply_page(&rows, q.page, q.page_size);
    Ok(Json(Page { rows: page, total }))
}

/// GET /runs/{id}/wga
pub async fn wga(
    State(state): State<SharedState>,
    Path(run_id): Path<i64>,
) -> ApiResult<Json<WgaData>> {
    let (project_id, query_ids, status) = jobs::run_meta(&state, run_id)?;
    if status != "succeeded" {
        return Err(ApiError::BadRequest(
            "This run has not finished yet. Please wait for it to complete.".into(),
        ));
    }
    let reference: Value = jobs::load_reference_json(&state, project_id, run_id)?;
    let lengths: Vec<(String, u64)> = serde_json::from_value(reference["lengths"].clone())
        .map_err(|_| ApiError::Internal("The reference metadata is unreadable.".into()))?;
    let genes: Vec<straincompass_types::WgaGene> =
        serde_json::from_value(reference["genes"].clone())
            .map_err(|_| ApiError::Internal("The reference metadata is unreadable.".into()))?;
    let mut queries = Vec::new();
    for qid in &query_ids {
        let res = jobs::load_query_result(&state, project_id, run_id, *qid)?;
        // Calls per gene, aligned with the reference gene list (keyed by
        // locus tag; genes missing from the coverage table count as absent).
        let by_tag: std::collections::HashMap<&str, &straincompass_types::GeneCoverageRow> = res
            .genes_coverage
            .iter()
            .map(|r| (r.locus_tag.as_str(), r))
            .collect();
        let mut calls = Vec::with_capacity(genes.len());
        let mut cov_pcts = Vec::with_capacity(genes.len());
        let mut identities = Vec::with_capacity(genes.len());
        for g in &genes {
            match by_tag.get(g.locus_tag.as_str()) {
                Some(r) => {
                    calls.push(r.call);
                    cov_pcts.push(r.cov_pct);
                    identities.push(r.best_identity);
                }
                None => {
                    calls.push(Call::Absent);
                    cov_pcts.push(0.0);
                    identities.push(0.0);
                }
            }
        }
        queries.push(WgaQuery {
            query_id: *qid,
            query_name: res.query_name.clone(),
            blocks: res.blocks,
            calls,
            cov_pcts,
            identities,
        });
    }
    Ok(Json(WgaData {
        reference: lengths,
        genes,
        queries,
    }))
}

/// GET /runs/{id}/gene/{locus} : preview stats + full MSA rows.
pub async fn gene_detail(
    State(state): State<SharedState>,
    Path((run_id, locus)): Path<(i64, String)>,
) -> ApiResult<Json<GeneDetail>> {
    let (project_id, _, _) = jobs::run_meta(&state, run_id)?;
    let (ref_fa, ref_gff, params, sources) = jobs::msa_sources(&state, project_id, run_id)?;
    let borrowed: Vec<_> = sources.iter().map(|s| s.borrow()).collect();
    let detail =
        straincompass_engine::pipeline::gene_detail(&ref_fa, &ref_gff, &params, &locus, &borrowed)?;
    Ok(Json(detail))
}

/// GET /runs/{id}/alignment : whole-genome alignment viewer data
/// (blocks + per-base variant events per query, events keyed by seqid).
///
/// Runs written recently carry variants.json already (run_comparison
/// precomputes the events), so this endpoint is normally a set of file
/// reads. Old runs backfill here: each query runs on its own blocking
/// task (bounded by cpu_slots) and the reference fasta is parsed once
/// for the whole run, instead of one serial loop that re-parsed the
/// reference under a global lock.
///
/// The answer is columnar (parallel arrays per event kind): the
/// object form was ~4 MB per query for divergent runs, and the
/// browser's JSON parser blocked the page for seconds on it.
pub async fn alignment(
    State(state): State<SharedState>,
    Path(run_id): Path<i64>,
) -> ApiResult<Json<straincompass_types::AlignmentDataColumnar>> {
    let (project_id, query_ids, status) = jobs::run_meta(&state, run_id)?;
    if status != "succeeded" {
        return Err(ApiError::BadRequest(
            "This run has not finished yet. Please wait for it to complete.".into(),
        ));
    }
    let reference: Value = jobs::load_reference_json(&state, project_id, run_id)?;
    let lengths: Vec<(String, u64)> = serde_json::from_value(reference["lengths"].clone())
        .map_err(|_| ApiError::Internal("The reference metadata is unreadable.".into()))?;

    // The parse fails with a friendly engine error when ref.fa is gone:
    // map it through like the job runner does.
    let ref_path = state
        .project_dir(project_id)
        .join("reference")
        .join("ref.fa");
    let ref_records = tokio::task::spawn_blocking(move || {
        straincompass_engine::fasta::parse_fasta(&ref_path).map_err(ApiError::from)
    })
    .await
    .map_err(|e| ApiError::Internal(format!("The alignment task crashed. ({e})")))??;
    let ref_records = std::sync::Arc::new(ref_records);

    let mut handles = Vec::with_capacity(query_ids.len());
    for qid in query_ids {
        let state = state.clone();
        let ref_records = std::sync::Arc::clone(&ref_records);
        let sem = state.cpu_slots.clone();
        handles.push(async move {
            let _permit = sem.acquire_owned().await;
            tokio::task::spawn_blocking(
                move || -> ApiResult<straincompass_types::AlignmentQueryColumnar> {
                    let (query_name, blocks) =
                        jobs::load_query_meta(&state, project_id, run_id, qid)?;
                    let events =
                        jobs::load_variants(&state, project_id, run_id, qid, &ref_records)?;
                    Ok(straincompass_types::AlignmentQueryColumnar {
                        query_id: qid,
                        query_name,
                        blocks,
                        events: events
                            .into_iter()
                            .map(|(seqid, ev)| (seqid, ev.into()))
                            .collect(),
                    })
                },
            )
            .await
            .map_err(|e| ApiError::Internal(format!("The alignment task crashed. ({e})")))?
        });
    }
    // Await in order: tasks were spawned eagerly so they run in
    // parallel; a later query's error does not hide an earlier one's.
    let mut queries = Vec::with_capacity(handles.len());
    for h in handles {
        queries.push(h.await?);
    }

    Ok(Json(straincompass_types::AlignmentDataColumnar {
        reference: lengths,
        queries,
    }))
}

#[derive(serde::Deserialize)]
pub struct RefseqQuery {
    pub seqid: String,
    pub start: u64,
    pub end: u64,
}

/// GET /runs/{id}/refseq?seqid&start&end : reference bases of a window,
/// for the alignment viewer's letters mode. Capped to keep responses
/// small; the frontend fetches only the visible window.
pub async fn refseq(
    State(state): State<SharedState>,
    Path(run_id): Path<i64>,
    Query(q): Query<RefseqQuery>,
) -> ApiResult<Json<Value>> {
    let (project_id, _, status) = jobs::run_meta(&state, run_id)?;
    if status != "succeeded" {
        return Err(ApiError::BadRequest(
            "This run has not finished yet. Please wait for it to complete.".into(),
        ));
    }
    const MAX_LEN: u64 = 8192;
    let start = q.start.max(1);
    let end = q.end.clamp(start, start + MAX_LEN - 1);
    let ref_fa = state
        .project_dir(project_id)
        .join("reference")
        .join("ref.fa");
    let records = straincompass_engine::fasta::parse_fasta(&ref_fa)?;
    let rec = records
        .iter()
        .find(|r| r.id == q.seqid)
        .ok_or_else(|| ApiError::NotFound("This reference sequence does not exist.".into()))?;
    let end = end.min(rec.seq.len() as u64);
    let seq = String::from_utf8_lossy(&straincompass_engine::fasta::subseq(rec, start, end, false))
        .into_owned();
    Ok(Json(serde_json::json!({
        "seqid": q.seqid,
        "start": start,
        "end": end,
        "seq": seq,
    })))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::state::AppState;
    use std::sync::{Arc, Mutex};
    use straincompass_types::RunParams;

    /// A tiny finished run on disk: one reference (chr1, 8 bp) and one
    /// query whose delta encodes an insertion, two deletions and a SNP
    /// (the same fixture the engine variant tests use).
    fn seeded_state() -> (SharedState, std::path::PathBuf) {
        let dir =
            std::env::temp_dir().join(format!("straincompass-api-test-{}", uuid::Uuid::new_v4()));
        let run_dir = dir.join("projects/1/runs/1");
        let qdir = run_dir.join("queries/10");
        std::fs::create_dir_all(qdir.join("work")).unwrap();
        std::fs::create_dir_all(dir.join("projects/1/reference")).unwrap();

        std::fs::write(dir.join("projects/1/reference/ref.fa"), ">chr1\nACGTACGT\n").unwrap();
        std::fs::write(
            run_dir.join("reference.json"),
            serde_json::to_vec(&serde_json::json!({
                "lengths": [["chr1", 8u64]],
                "genes": [{
                    "locus_tag": "gene1", "symbol": "gA", "biotype": "CDS",
                    "seqid": "chr1", "start": 1u64, "end": 8u64,
                    "strand": 1i8, "product": "test protein"
                }]
            }))
            .unwrap(),
        )
        .unwrap();
        std::fs::write(
            qdir.join("result.json"),
            serde_json::to_vec(&serde_json::json!({
                "query_name": "q1",
                "genes_coverage": [],
                "unaligned_gaps": [],
                "panel": null,
                "blocks": [{
                    "ref_seqid": "chr1", "ref_start": 1u64, "ref_end": 8u64,
                    "qry_seqid": "q1", "qry_start": 1u64, "qry_end": 6u64,
                    "qry_rev": false, "identity": 50.0f64
                }],
                "ref_lengths": [["chr1", 8u64]]
            }))
            .unwrap(),
        )
        .unwrap();
        std::fs::write(qdir.join("query.fa"), ">q1\nGGAAGT\n").unwrap();
        std::fs::write(
            qdir.join("work/cmp.delta"),
            ">chr1 q1\n1 8 1 6 1 0 0\n-1\n1\n1\n2\n0\n",
        )
        .unwrap();

        let conn = rusqlite::Connection::open_in_memory().unwrap();
        crate::db::init_db(&conn).unwrap();
        conn.execute("INSERT INTO projects (id, name) VALUES (1, 'test')", [])
            .unwrap();
        conn.execute(
            "INSERT INTO files (id, project_id, role, display_name, stored_name, size)
             VALUES (10, 1, 'query', 'q1.fasta', 'q1.fa', 12)",
            [],
        )
        .unwrap();
        conn.execute(
            "INSERT INTO runs (id, project_id, status, params_json, query_ids)
             VALUES (1, 1, 'succeeded', ?1, '[10]')",
            rusqlite::params![serde_json::to_string(&RunParams::default()).unwrap()],
        )
        .unwrap();

        let state: SharedState = Arc::new(AppState {
            db: Mutex::new(conn),
            data_dir: dir.clone(),
            cpu_slots: Arc::new(tokio::sync::Semaphore::new(1)),
        });
        (state, dir)
    }

    #[tokio::test]
    async fn alignment_returns_blocks_and_computes_events_on_demand() {
        let (state, dir) = seeded_state();
        let res = alignment(State(state.clone()), Path(1)).await.unwrap().0;

        assert_eq!(res.reference, vec![("chr1".to_string(), 8u64)]);
        assert_eq!(res.queries.len(), 1);
        let q = &res.queries[0];
        assert_eq!(q.query_name, "q1");
        assert_eq!(q.blocks.len(), 1);
        assert_eq!(q.blocks[0].ref_seqid, "chr1");

        let ev = &q.events["chr1"];
        assert_eq!(ev.ins_pos.len(), 1);
        assert_eq!((ev.ins_pos[0], ev.ins_seq[0].as_str()), (0, "G"));
        assert_eq!(ev.del_pos.len(), 2);
        assert_eq!((ev.del_pos[0], ev.del_len[0]), (1, 2));
        assert_eq!((ev.del_pos[1], ev.del_len[1]), (4, 1));
        assert_eq!(ev.snp_pos.len(), 1);
        assert_eq!(
            (ev.snp_pos[0], ev.snp_ref[0], ev.snp_qry[0]),
            (6, b'C', b'A')
        );

        // the events were cached next to result.json
        let cache = dir.join("projects/1/runs/1/queries/10/variants.json");
        assert!(cache.is_file());
        // a second call serves the same data from the cache
        let res2 = alignment(State(state), Path(1)).await.unwrap().0;
        assert_eq!(res2.queries[0].events["chr1"].snp_pos.len(), 1);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[tokio::test]
    async fn refseq_returns_the_requested_window() {
        let (state, dir) = seeded_state();
        let res = refseq(
            State(state),
            Path(1),
            Query(RefseqQuery {
                seqid: "chr1".to_string(),
                start: 3,
                end: 6,
            }),
        )
        .await
        .unwrap()
        .0;
        assert_eq!(res["seqid"], "chr1");
        assert_eq!(res["start"], 3);
        assert_eq!(res["end"], 6);
        assert_eq!(res["seq"], "GTAC");
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// Rewrite the seeded result.json with gained rows, as a run computed
    /// after the feature landed would have.
    fn seed_gained(dir: &std::path::Path) {
        let p = dir.join("projects/1/runs/1/queries/10/result.json");
        let mut v: serde_json::Value = serde_json::from_slice(&std::fs::read(&p).unwrap()).unwrap();
        v["gained"] = serde_json::json!([
            {
                "qry_seqid": "ctg1", "start": 1001u64, "end": 4000u64, "length": 3000u64,
                "gc_pct": 61.5f64, "at_contig_end": false, "anchor": "between",
                "anchor_seqid": "chr1", "anchor_start": 1000u64, "anchor_end": 1001u64,
                "flanks_disagree": false, "left_gene": "G1", "right_gene": "G2",
                "n_orfs": 3u32, "n_orfs_complete": 2u32, "orfs": []
            },
            {
                "qry_seqid": "ctg2", "start": 1u64, "end": 9000u64, "length": 9000u64,
                "gc_pct": 48.0f64, "at_contig_end": true, "anchor": "unanchored",
                "anchor_seqid": "", "anchor_start": 0u64, "anchor_end": 0u64,
                "flanks_disagree": false, "left_gene": "", "right_gene": "",
                "n_orfs": null, "n_orfs_complete": null, "orfs": []
            }
        ]);
        std::fs::write(&p, serde_json::to_vec(&v).unwrap()).unwrap();
    }

    fn gained_query(call: Option<&str>) -> TableQuery {
        TableQuery {
            query_id: None,
            page: 0,
            page_size: 200,
            sort_by: None,
            sort_dir: None,
            search: None,
            call: call.map(|s| s.to_string()),
            cols: None,
        }
    }

    /// Rewrite the seeded result.json with one gained region on the query
    /// contig that query.fa actually carries ("q1", 6 bp), so the verify
    /// endpoint can find both the row and its sequence.
    fn seed_gained_on_real_contig(dir: &std::path::Path) {
        let p = dir.join("projects/1/runs/1/queries/10/result.json");
        let mut v: serde_json::Value = serde_json::from_slice(&std::fs::read(&p).unwrap()).unwrap();
        v["gained"] = serde_json::json!([{
            "qry_seqid": "q1", "start": 1u64, "end": 6u64, "length": 6u64,
            "gc_pct": 50.0f64, "at_contig_end": true, "anchor": "unanchored",
            "anchor_seqid": "", "anchor_start": 0u64, "anchor_end": 0u64,
            "flanks_disagree": false, "left_gene": "", "right_gene": "",
            "n_orfs": 1u32, "n_orfs_complete": 1u32, "orfs": []
        }]);
        std::fs::write(&p, serde_json::to_vec(&v).unwrap()).unwrap();
    }

    /// Stand-in blast binaries, as the engine tests do for prodigal: the
    /// makeblastdb does nothing, the blastn writes a fixed hit table to
    /// whatever path follows -out. STRAINCOMPASS_TOOLS_DIRS points
    /// discovery at them, so no test needs BLAST+ installed.
    fn stub_blast_tools(dir: &std::path::Path) {
        use std::os::unix::fs::PermissionsExt;
        let bin = dir.join("stubbin");
        std::fs::create_dir_all(&bin).unwrap();
        // discovery is all-or-nothing, so the tools the back-check never
        // runs still have to exist by name.
        for name in ["nucmer", "show-coords", "show-snps", "dnadiff"] {
            std::fs::write(bin.join(name), "#!/bin/sh\nexit 1\n").unwrap();
        }
        std::fs::write(bin.join("makeblastdb"), "#!/bin/sh\nexit 0\n").unwrap();
        std::fs::write(
            bin.join("blastn"),
            "#!/bin/sh\nout=\"\"\nprev=\"\"\nfor a in \"$@\"; do\n  if [ \"$prev\" = \"-out\" ]; then out=\"$a\"; fi\n  prev=\"$a\"\ndone\ncat > \"$out\" <<'HITS'\nchr1\t100\t200\t99.5\t101\t1\t101\t1e-30\t185\nHITS\n",
        )
        .unwrap();
        for f in [
            "nucmer",
            "show-coords",
            "show-snps",
            "dnadiff",
            "makeblastdb",
            "blastn",
        ] {
            std::fs::set_permissions(bin.join(f), std::fs::Permissions::from_mode(0o755)).unwrap();
        }
        std::env::set_var("STRAINCOMPASS_TOOLS_DIRS", &bin);
    }

    #[tokio::test]
    async fn gained_verify_blats_the_region_back() {
        let (state, dir) = seeded_state();
        seed_gained_on_real_contig(&dir);
        stub_blast_tools(&dir);

        let v = gained_verify(
            State(state.clone()),
            Path(1),
            Query(GainedVerifyQuery {
                query_id: None,
                seqid: "q1".into(),
                start: 1,
                end: 6,
            }),
        )
        .await
        .unwrap()
        .0;
        assert_eq!(v.hits.len(), 1);
        assert_eq!(v.hits[0].ref_seqid, "chr1");
        assert_eq!((v.hits[0].ref_start, v.hits[0].ref_end), (100, 200));

        // A region that is not one of the run's rows is refused, not
        // searched: the endpoint reads the query fasta.
        let err = gained_verify(
            State(state.clone()),
            Path(1),
            Query(GainedVerifyQuery {
                query_id: None,
                seqid: "q1".into(),
                start: 2,
                end: 5,
            }),
        )
        .await
        .unwrap_err();
        assert!(matches!(err, ApiError::BadRequest(_)), "got {err:?}");

        std::env::remove_var("STRAINCOMPASS_TOOLS_DIRS");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[tokio::test]
    async fn gained_verify_asks_an_older_run_to_be_rerun() {
        let (state, dir) = seeded_state();
        stub_blast_tools(&dir);
        let err = gained_verify(
            State(state),
            Path(1),
            Query(GainedVerifyQuery {
                query_id: None,
                seqid: "q1".into(),
                start: 1,
                end: 6,
            }),
        )
        .await
        .unwrap_err();
        match err {
            ApiError::BadRequest(m) => assert!(m.contains("run the comparison again"), "got {m}"),
            other => panic!("expected BadRequest, got {other:?}"),
        }
        std::env::remove_var("STRAINCOMPASS_TOOLS_DIRS");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[tokio::test]
    async fn gained_returns_rows_and_filters_by_anchor() {
        let (state, dir) = seeded_state();
        seed_gained(&dir);

        let page = gained(State(state.clone()), Path(1), Query(gained_query(None)))
            .await
            .unwrap()
            .0;
        assert_eq!(page.total, 2);
        assert_eq!(
            page.rows[0].qry_seqid, "ctg1",
            "default sort is contig then start"
        );
        assert_eq!(page.rows[0].n_orfs, Some(3));
        // An absent count stays absent: it must never arrive as 0.
        assert_eq!(page.rows[1].n_orfs, None);

        let page = gained(
            State(state.clone()),
            Path(1),
            Query(gained_query(Some("unanchored"))),
        )
        .await
        .unwrap()
        .0;
        assert_eq!(page.total, 1);
        assert_eq!(page.rows[0].qry_seqid, "ctg2");

        let page = gained(State(state), Path(1), Query(gained_query(Some("anchored"))))
            .await
            .unwrap()
            .0;
        assert_eq!(page.total, 1);
        assert_eq!(page.rows[0].anchor, GainedAnchor::Between);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[tokio::test]
    async fn gained_asks_an_older_run_to_be_rerun() {
        // The seeded result.json predates the feature: no "gained" key at
        // all. It must still load (the other tables keep working) and this
        // endpoint must explain itself rather than fail.
        let (state, dir) = seeded_state();
        let err = gained(State(state), Path(1), Query(gained_query(None)))
            .await
            .unwrap_err();
        match err {
            ApiError::BadRequest(m) => assert!(m.contains("run the comparison again"), "got {m}"),
            other => panic!("expected BadRequest, got {other:?}"),
        }
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[tokio::test]
    async fn refseq_rejects_unknown_sequence() {
        let (state, dir) = seeded_state();
        let err = refseq(
            State(state),
            Path(1),
            Query(RefseqQuery {
                seqid: "nope".to_string(),
                start: 1,
                end: 5,
            }),
        )
        .await
        .unwrap_err();
        assert!(matches!(err, ApiError::NotFound(_)));
        let _ = std::fs::remove_dir_all(&dir);
    }
}
