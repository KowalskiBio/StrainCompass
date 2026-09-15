//! Strict BLAST recheck of a gene panel against each query genome.

use crate::tools::ToolPaths;
use crate::{friendly, Result};
use bactiment_types::{Call, PanelRow};
use std::collections::HashMap;
use std::io::Read;
use std::path::Path;
use std::process::Command;

/// Run the panel recheck for one query. `workdir` receives the blast
/// database files; it must be writable and private to this comparison.
#[allow(clippy::too_many_arguments)]
pub fn panel_recheck(
    tools: &ToolPaths,
    panel_fasta: &Path,
    query_fasta: &Path,
    workdir: &Path,
    query_name: &str,
    blast_cov: f64,
    blast_pid: f64,
    blast_evalue: f64,
) -> Result<Vec<PanelRow>> {
    std::fs::create_dir_all(workdir)?;
    let db = workdir.join("panel_db");

    let make_db = Command::new(&tools.makeblastdb)
        .args(["-in"])
        .arg(query_fasta)
        .args(["-dbtype", "nucl", "-out"])
        .arg(&db)
        .output()
        .map_err(|e| crate::EngineError::ToolMissing(format!("makeblastdb: {e}")))?;
    if !make_db.status.success() {
        return Err(friendly(format!(
            "The gene panel search could not be prepared for {}. {}",
            query_name,
            String::from_utf8_lossy(&make_db.stderr).trim()
        )));
    }

    let out = workdir.join("hits.tsv");
    let blast = Command::new(&tools.blastn)
        .args(["-query"])
        .arg(panel_fasta)
        .args(["-db"])
        .arg(&db)
        .args([
            "-outfmt",
            "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue qlen slen",
        ])
        .arg("-evalue")
        .arg(format!("{}", blast_evalue.max(1e-300)))
        .arg("-out")
        .arg(&out)
        .output()
        .map_err(|e| crate::EngineError::ToolMissing(format!("blastn: {e}")))?;
    if !blast.status.success() {
        return Err(friendly(format!(
            "The gene panel search failed for {}. {}",
            query_name,
            String::from_utf8_lossy(&blast.stderr).trim()
        )));
    }

    // Parse hits: per panel gene, merge query intervals for coverage,
    // identity weighted by aligned length.
    struct HitAcc {
        qlen: u64,
        ivs: Vec<(u64, u64)>,
        weighted_pid: f64,
        aligned: u64,
        best_evalue: f64,
    }
    let mut hits: HashMap<String, HitAcc> = HashMap::new();
    let mut text = String::new();
    std::fs::File::open(&out)?.read_to_string(&mut text)?;
    for line in text.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 13 {
            continue;
        }
        let gene_id = f[0].to_string();
        let pident: f64 = f[2].parse().unwrap_or(0.0);
        let length: u64 = f[3].parse().unwrap_or(0);
        let qstart: u64 = f[6].parse().unwrap_or(1);
        let qend: u64 = f[7].parse().unwrap_or(1);
        let evalue: f64 = f[10].parse().unwrap_or(f64::INFINITY);
        let qlen: u64 = f[11].parse().unwrap_or(0);
        let acc = hits.entry(gene_id).or_insert(HitAcc {
            qlen,
            ivs: Vec::new(),
            weighted_pid: 0.0,
            aligned: 0,
            best_evalue: f64::INFINITY,
        });
        acc.ivs.push((qstart.min(qend), qend.max(qstart)));
        acc.weighted_pid += pident * length as f64;
        acc.aligned += length;
        acc.best_evalue = acc.best_evalue.min(evalue);
    }

    let mut rows = Vec::new();
    let mut ids: Vec<String> = hits.keys().cloned().collect();
    ids.sort();
    for id in ids {
        let h = &hits[&id];
        let mut ivs = h.ivs.clone();
        ivs.sort();
        let mut merged: Vec<(u64, u64)> = Vec::new();
        for (s, e) in ivs {
            match merged.last_mut() {
                Some(last) if s <= last.1 + 1 => last.1 = last.1.max(e),
                _ => merged.push((s, e)),
            }
        }
        let cov_bp: u64 = merged.iter().map(|(s, e)| e - s + 1).sum();
        let cov_pct = if h.qlen > 0 {
            100.0 * cov_bp.min(h.qlen) as f64 / h.qlen as f64
        } else {
            0.0
        };
        let identity = if h.aligned > 0 {
            h.weighted_pid / h.aligned as f64
        } else {
            0.0
        };
        let call = if cov_pct >= blast_cov && identity >= blast_pid {
            Call::Present
        } else {
            Call::Absent
        };
        rows.push(PanelRow {
            gene_id: id,
            qlen: h.qlen,
            cov_pct,
            identity,
            best_evalue: format_evalue(h.best_evalue),
            call,
        });
    }
    // Panel genes without any hit: read the panel fasta for the id list.
    let panel_ids = panel_gene_ids(panel_fasta)?;
    for pid in panel_ids {
        if !hits.contains_key(&pid) {
            rows.push(PanelRow {
                gene_id: pid,
                qlen: 0,
                cov_pct: 0.0,
                identity: 0.0,
                best_evalue: "-".into(),
                call: Call::Absent,
            });
        }
    }
    rows.sort_by(|a, b| a.gene_id.cmp(&b.gene_id));
    Ok(rows)
}

fn format_evalue(e: f64) -> String {
    if e == f64::INFINITY {
        "-".into()
    } else if e == 0.0 {
        "0.0".into()
    } else if e >= 0.001 {
        format!("{e:.2}")
    } else {
        format!("{e:.1e}")
    }
}

fn panel_gene_ids(path: &Path) -> Result<Vec<String>> {
    let recs = crate::fasta::parse_fasta(path)?;
    Ok(recs.into_iter().map(|r| r.id).collect())
}
