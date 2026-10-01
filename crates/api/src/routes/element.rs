//! Panel genes beyond the call: the panel across all queries of a run,
//! where a gene sits in one query (plasmid, chromosomal insertion,
//! mobile-element neighbours), and how much of the element carrying it
//! every other query holds.

use crate::error::{ApiError, ApiResult};
use crate::jobs;
use crate::routes::results::{apply_page, gained_examine, matches_search, resolve_query_id};
use crate::state::SharedState;
use axum::extract::{Path, Query, State};
use axum::Json;
use straincompass_engine::delta::DeltaFile;
use straincompass_engine::element::{self, ElementTarget};
use straincompass_engine::fasta;
use straincompass_types::{
    Call, ElementReport, Page, PanelContext, PanelMatrixRow, PanelRow, TableQuery,
};

fn succeeded_queries(state: &SharedState, run_id: i64) -> ApiResult<(i64, Vec<i64>)> {
    let (project_id, query_ids, status) = jobs::run_meta(state, run_id)?;
    if status != "succeeded" {
        return Err(ApiError::BadRequest(
            "This run has not finished yet. Please wait for it to complete.".into(),
        ));
    }
    Ok((project_id, query_ids))
}

fn no_panel() -> ApiError {
    ApiError::BadRequest("This run has no gene panel results (no panel was provided).".into())
}

/// GET /runs/{id}/panel_matrix : every panel gene across all queries.
pub async fn panel_matrix(
    State(state): State<SharedState>,
    Path(run_id): Path<i64>,
    Query(q): Query<TableQuery>,
) -> ApiResult<Json<Page<PanelMatrixRow>>> {
    let (project_id, query_ids) = succeeded_queries(&state, run_id)?;
    let mut panels: Vec<Vec<PanelRow>> = Vec::with_capacity(query_ids.len());
    for qid in &query_ids {
        let res = jobs::load_query_result(&state, project_id, run_id, *qid)?;
        panels.push(res.panel.ok_or_else(no_panel)?);
    }
    let mut rows: Vec<PanelMatrixRow> = Vec::new();
    if let Some(first) = panels.first() {
        for g in first {
            let mut row = PanelMatrixRow {
                gene_id: g.gene_id.clone(),
                qlen: g.qlen,
                ..Default::default()
            };
            for panel in &panels {
                // every query was searched with the same panel, but look
                // genes up by id rather than trusting the order
                let hit = panel.iter().find(|r| r.gene_id == g.gene_id);
                row.calls.push(hit.map(|h| h.call).unwrap_or(Call::Absent));
                row.cov_pcts.push(hit.map(|h| h.cov_pct).unwrap_or(0.0));
                row.identities.push(hit.map(|h| h.identity).unwrap_or(0.0));
                row.loci
                    .push(hit.map(|h| h.qry_locus.clone()).unwrap_or_default());
            }
            rows.push(row);
        }
    }
    rows.retain(|r| matches_search(&r.gene_id, &q.search));
    if q.call.as_deref() == Some("not_present") {
        rows.retain(|r| r.calls.iter().any(|c| *c != Call::Present));
    }
    let asc = q.sort_dir.as_deref() != Some("desc");
    rows.sort_by(|a, b| {
        let ord = match q.sort_by.as_deref() {
            Some("qlen") => a.qlen.cmp(&b.qlen),
            _ => a.gene_id.to_lowercase().cmp(&b.gene_id.to_lowercase()),
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

#[derive(serde::Deserialize)]
pub struct PanelContextQuery {
    pub query_id: Option<i64>,
    pub gene_id: String,
}

/// GET /runs/{id}/panel_context?query_id&gene_id : where a panel gene
/// sits in one query.
pub async fn panel_context(
    State(state): State<SharedState>,
    Path(run_id): Path<i64>,
    Query(q): Query<PanelContextQuery>,
) -> ApiResult<Json<PanelContext>> {
    let (project_id, run_id, qid) = resolve_query_id(&state, run_id, q.query_id)?;
    let res = jobs::load_query_result(&state, project_id, run_id, qid)?;
    let panel = res.panel.ok_or_else(no_panel)?;
    let row = panel
        .into_iter()
        .find(|r| r.gene_id == q.gene_id)
        .ok_or_else(|| ApiError::NotFound(format!("{} is not in this run's panel.", q.gene_id)))?;
    let qdir = state
        .run_dir(project_id, run_id)
        .join("queries")
        .join(qid.to_string());
    let genes_note = match &res.gained_orfs {
        straincompass_types::GainedOrfStatus::Unavailable(reason) => reason.clone(),
        _ => String::new(),
    };
    let mut gained = res.gained.unwrap_or_default();
    crate::routes::nblast::merge_ncbi_names(&qdir, &mut gained);
    let query_name = res.query_name;
    // the reference's own genes, to list neighbours in shared stretches;
    // a run without them still answers, with predicted genes only
    let ref_genes: Vec<straincompass_types::WgaGene> =
        jobs::load_reference_json(&state, project_id, run_id)
            .ok()
            .and_then(|v| serde_json::from_value(v["genes"].clone()).ok())
            .unwrap_or_default();

    let ctx = tokio::task::spawn_blocking(move || -> ApiResult<PanelContext> {
        let qry_fa = qdir.join("query.fa");
        let delta_path = qdir.join("work").join("cmp.delta");
        if !qry_fa.is_file() || !delta_path.is_file() {
            return Err(ApiError::NotFound(
                "This run's alignment files are no longer on the server, so the gene's surroundings cannot be shown.".into(),
            ));
        }
        let records = fasta::parse_fasta(&qry_fa)?;
        let delta = DeltaFile::parse(&delta_path)?;
        let stats = element::contig_stats(&records, &delta);
        let mut ctx = PanelContext {
            query_id: qid,
            query_name,
            gene_id: row.gene_id.clone(),
            call: row.call,
            cov_pct: row.cov_pct,
            identity: row.identity,
            n_contigs: stats.len(),
            genome_bp: stats.iter().map(|s| s.length).sum(),
            window: element::CONTEXT_WINDOW,
            ..Default::default()
        };
        let Some((contig, start, end, strand)) = element::parse_locus(&row.qry_locus) else {
            ctx.verdict = "not_found".into();
            ctx.verdict_text = if row.call == Call::Absent {
                format!("{} was not found in this genome.", row.gene_id)
            } else {
                "This run did not keep the hit's position; run the comparison again to see where the gene sits.".into()
            };
            return Ok(ctx);
        };
        let stat = stats
            .iter()
            .find(|s| s.seqid == contig)
            .cloned()
            .unwrap_or_default();
        let region = gained
            .iter()
            .find(|r| r.qry_seqid == contig && r.start <= end && r.end >= start)
            .cloned();
        let mut genes = element::context_genes(
            &gained,
            &contig,
            (start, end),
            element::CONTEXT_WINDOW,
            &row.gene_id,
        );
        genes.extend(element::annotation_context_genes(
            &delta.alignments,
            &ref_genes,
            &contig,
            (start, end),
            element::CONTEXT_WINDOW,
        ));
        genes.sort_by_key(|g| g.start);
        let mut markers: Vec<String> = Vec::new();
        for g in genes.iter().filter(|g| g.mobile && !g.is_hit) {
            if !markers.contains(&g.label) {
                markers.push(g.label.clone());
            }
        }
        let (verdict, text) = element::verdict(&stat, region.as_ref(), &markers);
        ctx.contig = Some(stat);
        ctx.hit_start = start;
        ctx.hit_end = end;
        ctx.hit_strand = strand;
        ctx.region = region;
        ctx.genes = genes;
        ctx.verdict = verdict;
        ctx.verdict_text = text;
        ctx.mobile_markers = markers;
        ctx.genes_note = genes_note;
        Ok(ctx)
    })
    .await
    .map_err(|e| ApiError::Internal(format!("The gene's surroundings could not be read. ({e})")))??;
    Ok(Json(ctx))
}

#[derive(serde::Deserialize)]
pub struct PanelElementBody {
    pub gene_id: String,
    pub source_query_id: i64,
    /// A complete record (e.g. a plasmid) to compare instead of the
    /// contig that carries the gene in the source query.
    #[serde(default)]
    pub accession: Option<String>,
}

/// POST /runs/{id}/panel_element : how much of the element carrying a
/// panel gene every query of the run holds.
pub async fn panel_element(
    State(state): State<SharedState>,
    Path(run_id): Path<i64>,
    Json(body): Json<PanelElementBody>,
) -> ApiResult<Json<ElementReport>> {
    let (project_id, query_ids) = succeeded_queries(&state, run_id)?;
    if !query_ids.contains(&body.source_query_id) {
        return Err(ApiError::BadRequest(
            "This query is not part of the run.".into(),
        ));
    }
    let run_dir = state.run_dir(project_id, run_id);
    let qdir = |qid: i64| run_dir.join("queries").join(qid.to_string());

    // names and calls per query
    let mut meta: Vec<(i64, String, Option<Call>, Option<String>)> = Vec::new();
    for qid in &query_ids {
        let res = jobs::load_query_result(&state, project_id, run_id, *qid)?;
        let hit = res
            .panel
            .as_ref()
            .and_then(|p| p.iter().find(|r| r.gene_id == body.gene_id));
        meta.push((
            *qid,
            res.query_name,
            hit.map(|h| h.call),
            hit.map(|h| h.qry_locus.clone()),
        ));
    }

    let accession = body
        .accession
        .as_deref()
        .map(str::trim)
        .filter(|a| !a.is_empty())
        .map(str::to_string);
    let (kind, name, title, fetched) = match &accession {
        Some(acc) => {
            let key = {
                let conn = state.db.lock().unwrap();
                crate::db::get_setting(&conn, "ncbi_api_key")?
            };
            let (title, seq) = crate::routes::ncbi::fetch_record(acc, key.as_deref()).await?;
            ("accession", acc.clone(), title, Some(seq))
        }
        None => {
            let locus = meta
                .iter()
                .find(|m| m.0 == body.source_query_id)
                .and_then(|m| m.3.clone())
                .unwrap_or_default();
            let contig = element::parse_locus(&locus)
                .map(|l| l.0)
                .or_else(|| (!locus.is_empty()).then(|| locus.clone()))
                .ok_or_else(|| {
                    ApiError::BadRequest(format!(
                        "{} was not found in the chosen strain, so there is no contig to compare. Choose a strain where it is present.",
                        body.gene_id
                    ))
                })?;
            ("contig", contig, String::new(), None)
        }
    };

    let src_fa = qdir(body.source_query_id).join("query.fa");
    let targets: Vec<ElementTarget> = query_ids
        .iter()
        .map(|q| ElementTarget {
            query_id: *q,
            fasta: qdir(*q).join("query.fa"),
            db: Some(qdir(*q).join("work").join("panel").join("panel_db")),
        })
        .collect();
    if targets.iter().any(|t| !t.fasta.is_file()) || (fetched.is_none() && !src_fa.is_file()) {
        return Err(ApiError::NotFound(
            "This run's genome files are no longer on the server, so the element cannot be compared.".into(),
        ));
    }

    let element_name = name.clone();
    let (element_len, sizes, covs) = gained_examine(
        &state,
        project_id,
        run_id,
        body.source_query_id,
        move |tools, work| {
            let seq = match fetched {
                Some(s) => s,
                None => fasta::parse_fasta(&src_fa)?
                    .into_iter()
                    .find(|r| r.id == element_name)
                    .map(|r| r.seq)
                    .ok_or_else(|| {
                        straincompass_engine::friendly(format!(
                            "The contig {element_name} is not in the strain's genome file."
                        ))
                    })?,
            };
            let sizes: Vec<u64> = targets
                .iter()
                .map(|t| {
                    fasta::parse_fasta(&t.fasta)
                        .map(|rs| rs.iter().map(|r| r.seq.len() as u64).sum())
                })
                .collect::<straincompass_engine::Result<_>>()?;
            let covs = element::element_presence(tools, &element_name, &seq, &targets, &work)?;
            Ok((seq.len() as u64, sizes, covs))
        },
    )
    .await?;

    let hits = covs
        .into_iter()
        .zip(sizes)
        .map(|((qid, cov), genome_bp)| {
            let m = meta.iter().find(|m| m.0 == qid);
            element::element_hit(
                qid,
                m.map(|m| m.1.clone()).unwrap_or_default(),
                m.and_then(|m| m.2),
                cov,
                element_len,
                genome_bp,
            )
        })
        .collect();
    Ok(Json(ElementReport {
        gene_id: body.gene_id,
        source_query_id: body.source_query_id,
        element_kind: kind.into(),
        element_name: name,
        element_title: title,
        element_len,
        hits,
    }))
}
