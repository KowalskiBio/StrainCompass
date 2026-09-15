//! Result table endpoints with server side sort/filter/pagination,
//! the WGA viewer data and the gene alignment (MSA) endpoint.

use crate::error::{ApiError, ApiResult};
use crate::jobs;
use crate::state::SharedState;
use axum::extract::{Path, Query, State};
use axum::Json;
use bactiment_types::{
    Call, GapRow, GeneCoverageRow, GeneDetail, MatrixRow, Page, PanelRow, TableQuery,
    WgaData, WgaQuery,
};
use serde_json::Value;

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
        None => *query_ids.first().ok_or_else(|| {
            ApiError::BadRequest("This run has no queries.".into())
        })?,
    };
    Ok((project_id, run_id, qid))
}

fn matches_search(row_text: &str, search: &Option<String>) -> bool {
    match search {
        Some(s) if !s.trim().is_empty() => row_text.to_lowercase().contains(&s.trim().to_lowercase()),
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

fn sort_rows(rows: &mut Vec<GeneCoverageRow>, by: &Option<String>, dir: &Option<String>) {
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
            "cov_pct" => a.cov_pct.partial_cmp(&b.cov_pct).unwrap_or(std::cmp::Ordering::Equal),
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
    rows.retain(|r| {
        matches_search(&row_str(r), &q.search) && call_filter_ok(r.call, &q.call)
    });
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

/// GET /runs/{id}/panel_recheck
pub async fn panel_recheck(
    State(state): State<SharedState>,
    Path(run_id): Path<i64>,
    Query(q): Query<TableQuery>,
) -> ApiResult<Json<Page<PanelRow>>> {
    let (project_id, run_id, qid) = resolve_query_id(&state, run_id, q.query_id)?;
    let res = jobs::load_query_result(&state, project_id, run_id, qid)?;
    let mut rows = res.panel.ok_or_else(|| {
        ApiError::BadRequest(
            "This run has no gene panel results (no panel was provided).".into(),
        )
    })?;
    rows.retain(|r| {
        matches_search(&r.gene_id, &q.search) && call_filter_ok(r.call, &q.call)
    });
    let asc = q.sort_dir.as_deref() != Some("desc");
    let by = q.sort_by.clone().unwrap_or_default();
    rows.sort_by(|a, b| {
        let ord = match by.as_str() {
            "gene_id" => a.gene_id.cmp(&b.gene_id),
            "qlen" => a.qlen.cmp(&b.qlen),
            "cov_pct" => a.cov_pct.partial_cmp(&b.cov_pct).unwrap_or(std::cmp::Ordering::Equal),
            "identity" => a.identity.partial_cmp(&b.identity).unwrap_or(std::cmp::Ordering::Equal),
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
    rows.retain(|r| matches_search(&format!("{} {} {}", r.locus_tag, r.symbol, r.seqid), &q.search));
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
pub async fn wga(State(state): State<SharedState>, Path(run_id): Path<i64>) -> ApiResult<Json<WgaData>> {
    let (project_id, query_ids, status) = jobs::run_meta(&state, run_id)?;
    if status != "succeeded" {
        return Err(ApiError::BadRequest(
            "This run has not finished yet. Please wait for it to complete.".into(),
        ));
    }
    let reference: Value = jobs::load_reference_json(&state, project_id, run_id)?;
    let lengths: Vec<(String, u64)> = serde_json::from_value(reference["lengths"].clone())
        .map_err(|_| ApiError::Internal("The reference metadata is unreadable.".into()))?;
    let genes: Vec<bactiment_types::WgaGene> = serde_json::from_value(reference["genes"].clone())
        .map_err(|_| ApiError::Internal("The reference metadata is unreadable.".into()))?;
    let mut queries = Vec::new();
    for qid in &query_ids {
        let res = jobs::load_query_result(&state, project_id, run_id, *qid)?;
        queries.push(WgaQuery {
            query_id: *qid,
            query_name: res.query_name.clone(),
            blocks: res.blocks,
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
    let detail = bactiment_engine::pipeline::gene_detail(
        &ref_fa,
        &ref_gff,
        &params,
        &locus,
        &borrowed,
    )?;
    Ok(Json(detail))
}

