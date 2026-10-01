//! The resistance and virulence screen results: per query, across all
//! queries of a run, and whether the screen ran for each query.

use crate::error::{ApiError, ApiResult};
use crate::jobs;
use crate::routes::results::{apply_page, matches_search, resolve_query_id};
use crate::state::SharedState;
use axum::extract::{Path, Query, State};
use axum::Json;
use std::collections::BTreeMap;
use straincompass_types::{Page, ScreenHit, ScreenMatrixRow, ScreenStatus, TableQuery};

/// The `call` filter of the screen tables picks a kind of gene.
fn kind_ok(kind: &str, filter: &Option<String>) -> bool {
    match filter.as_deref() {
        Some(f) if !f.is_empty() && f != "not_present" => kind.eq_ignore_ascii_case(f),
        _ => true,
    }
}

fn hit_text(h: &ScreenHit) -> String {
    format!(
        "{} {} {} {} {} {}",
        h.gene, h.product, h.kind, h.category, h.class, h.contig
    )
}

/// GET /runs/{id}/screen?query_id : one query's hits.
pub async fn screen(
    State(state): State<SharedState>,
    Path(run_id): Path<i64>,
    Query(q): Query<TableQuery>,
) -> ApiResult<Json<Page<ScreenHit>>> {
    let (project_id, run_id, qid) = resolve_query_id(&state, run_id, q.query_id)?;
    let res = jobs::load_query_result(&state, project_id, run_id, qid)?;
    let mut rows = res.screen.unwrap_or_default();
    rows.retain(|h| matches_search(&hit_text(h), &q.search) && kind_ok(&h.kind, &q.call));
    let asc = q.sort_dir.as_deref() != Some("desc");
    rows.sort_by(|a, b| {
        let ord = match q.sort_by.as_deref() {
            Some("gene") => a.gene.to_lowercase().cmp(&b.gene.to_lowercase()),
            Some("kind") => {
                (a.kind.as_str(), a.category.as_str()).cmp(&(b.kind.as_str(), b.category.as_str()))
            }
            Some("class") => a.class.cmp(&b.class),
            Some("source") => a.source.cmp(&b.source),
            Some("identity") => a
                .identity
                .partial_cmp(&b.identity)
                .unwrap_or(std::cmp::Ordering::Equal),
            Some("coverage") => a
                .coverage
                .partial_cmp(&b.coverage)
                .unwrap_or(std::cmp::Ordering::Equal),
            _ => (a.kind.as_str(), a.category.as_str(), a.gene.to_lowercase()).cmp(&(
                b.kind.as_str(),
                b.category.as_str(),
                b.gene.to_lowercase(),
            )),
        };
        if asc {
            ord
        } else {
            ord.reverse()
        }
    });
    let total = rows.len() as u64;
    Ok(Json(Page {
        rows: apply_page(&rows, q.page, q.page_size),
        total,
    }))
}

#[derive(serde::Serialize)]
pub struct QueryScreenStatus {
    pub query_id: i64,
    pub status: ScreenStatus,
}

/// GET /runs/{id}/screen_status : whether the screen ran, per query.
pub async fn screen_status(
    State(state): State<SharedState>,
    Path(run_id): Path<i64>,
) -> ApiResult<Json<Vec<QueryScreenStatus>>> {
    let (project_id, query_ids, status) = jobs::run_meta(&state, run_id)?;
    if status != "succeeded" {
        return Err(ApiError::BadRequest(
            "This run has not finished yet. Please wait for it to complete.".into(),
        ));
    }
    let mut out = Vec::new();
    for qid in query_ids {
        let res = jobs::load_query_result(&state, project_id, run_id, qid)?;
        out.push(QueryScreenStatus {
            query_id: qid,
            status: res.screen_status,
        });
    }
    Ok(Json(out))
}

/// GET /runs/{id}/screen_matrix : every screened gene across all queries.
pub async fn screen_matrix(
    State(state): State<SharedState>,
    Path(run_id): Path<i64>,
    Query(q): Query<TableQuery>,
) -> ApiResult<Json<Page<ScreenMatrixRow>>> {
    let (project_id, query_ids, status) = jobs::run_meta(&state, run_id)?;
    if status != "succeeded" {
        return Err(ApiError::BadRequest(
            "This run has not finished yet. Please wait for it to complete.".into(),
        ));
    }
    let n = query_ids.len();
    // one row per (source, gene): a gene found by both databases shows
    // twice, each with its own evidence
    let mut rows: BTreeMap<(String, String, String), ScreenMatrixRow> = BTreeMap::new();
    for (i, qid) in query_ids.iter().enumerate() {
        let res = jobs::load_query_result(&state, project_id, run_id, *qid)?;
        for h in res.screen.unwrap_or_default() {
            let key = (h.kind.clone(), h.source.clone(), h.gene.clone());
            let row = rows.entry(key).or_insert_with(|| ScreenMatrixRow {
                source: h.source.clone(),
                gene: h.gene.clone(),
                product: h.product.clone(),
                kind: h.kind.clone(),
                category: h.category.clone(),
                class: h.class.clone(),
                identities: vec![None; n],
            });
            // several copies in one genome: the best one stands for it
            let best = row.identities[i].map_or(h.identity, |v| v.max(h.identity));
            row.identities[i] = Some(best);
        }
    }
    let mut rows: Vec<ScreenMatrixRow> = rows.into_values().collect();
    rows.retain(|r| {
        matches_search(
            &format!(
                "{} {} {} {} {}",
                r.gene, r.product, r.kind, r.category, r.class
            ),
            &q.search,
        ) && kind_ok(&r.kind, &q.call)
    });
    if q.call.as_deref() == Some("not_present") {
        rows.retain(|r| r.identities.iter().any(|v| v.is_none()));
    }
    let total = rows.len() as u64;
    Ok(Json(Page {
        rows: apply_page(&rows, q.page, q.page_size),
        total,
    }))
}
