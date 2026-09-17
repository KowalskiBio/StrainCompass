//! TSV/CSV export endpoints. Exports respect the current view: sorting,
//! filters and search are applied exactly as in the table endpoints.

use crate::error::{ApiError, ApiResult};
use crate::jobs;
use crate::state::SharedState;
use axum::body::Body;
use axum::extract::{Path, Query, State};
use axum::http::header;
use axum::http::HeaderMap;
use axum::response::{IntoResponse, Response};
use serde::Deserialize;
use straincompass_types::{Call, GeneCoverageRow, TableQuery};

#[derive(Debug, Deserialize)]
pub struct ExportQuery {
    #[serde(flatten)]
    pub table: TableQuery,
    #[serde(default)]
    pub format: Option<String>,
}

fn attachment_headers(file_name: &str, csv: bool) -> HeaderMap {
    let mut h = HeaderMap::new();
    h.insert(
        header::CONTENT_DISPOSITION,
        format!("attachment; filename=\"{file_name}\"")
            .parse()
            .unwrap(),
    );
    h.insert(
        header::CONTENT_TYPE,
        if csv {
            "text/csv; charset=utf-8".parse().unwrap()
        } else {
            "text/tab-separated-values; charset=utf-8".parse().unwrap()
        },
    );
    h
}

fn csv_escape(s: &str) -> String {
    if s.contains(',') || s.contains('"') || s.contains('\n') {
        format!("\"{}\"", s.replace('"', "\"\""))
    } else {
        s.to_string()
    }
}

fn fmt_num(v: f64) -> String {
    format!("{v:.2}")
}

/// GET /runs/{id}/export/{table}?format=tsv|csv&cols=...
pub async fn export_table(
    State(state): State<SharedState>,
    Path((run_id, table)): Path<(i64, String)>,
    Query(eq): Query<ExportQuery>,
) -> ApiResult<Response> {
    let q = eq.table;
    let csv = eq.format.as_deref() == Some("csv");

    let (project_id, query_ids, status) = jobs::run_meta(&state, run_id)?;
    if status != "succeeded" {
        return Err(ApiError::BadRequest(
            "This run has not finished yet. Please wait for it to complete.".into(),
        ));
    }
    let qid = q
        .query_id
        .unwrap_or_else(|| query_ids.first().copied().unwrap_or(0));
    let sep = if csv { "," } else { "\t" };
    let (filename, body): (String, String) = match table.as_str() {
        "genes_coverage" => {
            let res = jobs::load_query_result(&state, project_id, run_id, qid)?;
            let mut rows = res.genes_coverage;
            rows.retain(|r| filter_row(r, &q));
            sort_coverage(&mut rows, &q);
            let cols = coverage_cols(&q.cols);
            let mut lines = vec![header_line(&cols, sep)];
            for r in &rows {
                lines.push(coverage_line(r, &cols, sep));
            }
            (
                format!("{}_genes_coverage.{}", res.query_name, if csv { "csv" } else { "tsv" }),
                lines.join("\n") + "\n",
            )
        }
        "unaligned_gaps" => {
            let res = jobs::load_query_result(&state, project_id, run_id, qid)?;
            let mut rows = res.unaligned_gaps;
            rows.retain(|r| {
                let text = format!("{} {}", r.seqid, r.genes.join(" "));
                match &q.search {
                    Some(s) if !s.trim().is_empty() => {
                        text.to_lowercase().contains(&s.trim().to_lowercase())
                    }
                    _ => true,
                }
            });
            let mut lines = vec![join(&[
                "seqid".into(),
                "start".into(),
                "end".into(),
                "length".into(),
                "n_genes".into(),
                "genes".into(),
            ], sep)];
            for r in &rows {
                lines.push(join(
                    &[
                        r.seqid.clone(),
                        r.start.to_string(),
                        r.end.to_string(),
                        r.length.to_string(),
                        r.n_genes.to_string(),
                        r.genes.join(";"),
                    ],
                    sep,
                ));
            }
            (
                format!("{}_unaligned_gaps.{}", res.query_name, if csv { "csv" } else { "tsv" }),
                lines.join("\n") + "\n",
            )
        }
        "panel_recheck" => {
            let res = jobs::load_query_result(&state, project_id, run_id, qid)?;
            let mut rows = res.panel.ok_or_else(|| {
                ApiError::BadRequest("This run has no gene panel results.".into())
            })?;
            rows.retain(|r| match &q.search {
                Some(s) if !s.trim().is_empty() => {
                    r.gene_id.to_lowercase().contains(&s.trim().to_lowercase())
                }
                _ => true,
            });
            let mut lines = vec![join(
                &[
                    "gene_id".into(),
                    "qlen".into(),
                    "cov_pct".into(),
                    "identity".into(),
                    "best_evalue".into(),
                    "call".into(),
                ],
                sep,
            )];
            for r in &rows {
                lines.push(join(
                    &[
                        r.gene_id.clone(),
                        r.qlen.to_string(),
                        fmt_num(r.cov_pct),
                        fmt_num(r.identity),
                        r.best_evalue.clone(),
                        r.call.as_str().to_string(),
                    ],
                    sep,
                ));
            }
            (
                format!("{}_panel_recheck.{}", res.query_name, if csv { "csv" } else { "tsv" }),
                lines.join("\n") + "\n",
            )
        }
        "matrix" => {
            let mut results = Vec::new();
            for qid in &query_ids {
                results.push(jobs::load_query_result(&state, project_id, run_id, *qid)?);
            }
            let mut header = vec![
                "locus_tag".to_string(),
                "symbol".to_string(),
                "biotype".to_string(),
                "seqid".to_string(),
                "start".to_string(),
                "end".to_string(),
            ];
            for res in &results {
                header.push(res.query_name.clone());
            }
            let mut lines = vec![join(&header, sep)];
            if let Some(first) = results.first() {
                for (i, g) in first.genes_coverage.iter().enumerate() {
                    let text = format!("{} {} {}", g.locus_tag, g.symbol, g.seqid);
                    let ok = match &q.search {
                        Some(s) if !s.trim().is_empty() => {
                            text.to_lowercase().contains(&s.trim().to_lowercase())
                        }
                        _ => true,
                    };
                    if !ok {
                        continue;
                    }
                    let mut fields = vec![
                        g.locus_tag.clone(),
                        g.symbol.clone(),
                        g.biotype.clone(),
                        g.seqid.clone(),
                        g.start.to_string(),
                        g.end.to_string(),
                    ];
                    for res in &results {
                        fields.push(res.genes_coverage[i].call.as_str().to_string());
                    }
                    lines.push(join(&fields, sep));
                }
            }
            (
                format!("presence_absence_matrix.{}", if csv { "csv" } else { "tsv" }),
                lines.join("\n") + "\n",
            )
        }
        _ => {
            return Err(ApiError::BadRequest(
                "Unknown table for export. Choose genes_coverage, unaligned_gaps, panel_recheck or matrix.".into(),
            ))
        }
    };
    let headers = attachment_headers(&filename, csv);
    Ok((headers, Body::from(body)).into_response())
}

fn join(fields: &[String], sep: &str) -> String {
    if sep == "," {
        fields
            .iter()
            .map(|f| csv_escape(f))
            .collect::<Vec<_>>()
            .join(",")
    } else {
        fields.join("\t")
    }
}

fn filter_row(r: &GeneCoverageRow, q: &TableQuery) -> bool {
    let text = format!("{} {} {}", r.locus_tag, r.symbol, r.seqid);
    let search_ok = match &q.search {
        Some(s) if !s.trim().is_empty() => text.to_lowercase().contains(&s.trim().to_lowercase()),
        _ => true,
    };
    let call_ok = match q.call.as_deref() {
        Some("present") => r.call == Call::Present,
        Some("partial") => r.call == Call::Partial,
        Some("absent") => r.call == Call::Absent,
        _ => true,
    };
    search_ok && call_ok
}

fn sort_coverage(rows: &mut [GeneCoverageRow], q: &TableQuery) {
    let asc = q.sort_dir.as_deref() != Some("desc");
    let by = q.sort_by.clone().unwrap_or_default();
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

const COVERAGE_ALL_COLS: &[(&str, &str)] = &[
    ("locus_tag", "Locus tag"),
    ("symbol", "Symbol"),
    ("protein_id", "Protein accession"),
    ("biotype", "Biotype"),
    ("seqid", "Sequence"),
    ("start", "Start"),
    ("end", "End"),
    ("length", "Length"),
    ("cov_bp", "Covered bases"),
    ("cov_pct", "Coverage %"),
    ("call", "Call"),
    ("best_identity", "Best identity %"),
    ("mismatches", "Mismatches"),
    ("indels", "Indels"),
];

fn coverage_cols(cols: &Option<String>) -> Vec<String> {
    match cols {
        Some(c) if !c.trim().is_empty() => c
            .split(',')
            .map(|s| s.trim().to_string())
            .filter(|s| COVERAGE_ALL_COLS.iter().any(|(k, _)| *k == *s))
            .collect(),
        _ => COVERAGE_ALL_COLS
            .iter()
            .map(|(k, _)| k.to_string())
            .collect(),
    }
}

fn header_line(cols: &[String], sep: &str) -> String {
    let labels: Vec<String> = cols
        .iter()
        .map(|c| {
            COVERAGE_ALL_COLS
                .iter()
                .find(|(k, _)| k == c)
                .map(|(_, l)| l.to_string())
                .unwrap_or_else(|| c.clone())
        })
        .collect();
    join(&labels, sep)
}

fn coverage_line(r: &GeneCoverageRow, cols: &[String], sep: &str) -> String {
    let fields: Vec<String> = cols
        .iter()
        .map(|c| match c.as_str() {
            "locus_tag" => r.locus_tag.clone(),
            "symbol" => r.symbol.clone(),
            "protein_id" => r.protein_id.clone(),
            "biotype" => r.biotype.clone(),
            "seqid" => r.seqid.clone(),
            "start" => r.start.to_string(),
            "end" => r.end.to_string(),
            "length" => r.length.to_string(),
            "cov_bp" => r.cov_bp.to_string(),
            "cov_pct" => fmt_num(r.cov_pct),
            "call" => r.call.as_str().to_string(),
            "best_identity" => fmt_num(r.best_identity),
            "mismatches" => r.mismatches.to_string(),
            "indels" => r.indels.to_string(),
            _ => String::new(),
        })
        .collect();
    join(&fields, sep)
}
