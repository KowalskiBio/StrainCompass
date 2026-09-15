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

    // Parse hits. R parity: the best hit per panel gene is the single
    // row with the highest bitscore; identity and coverage are that
    // row's pident and qcovs, not a merge over all HSPs.
    let out = workdir.join("hits.tsv");
    let blast = Command::new(&tools.blastn)
        .args(["-query"])
        .arg(panel_fasta)
        .args(["-db"])
        .arg(&db)
        .args([
            "-outfmt",
            "6 qseqid sseqid pident length qlen qcovs evalue bitscore",
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

    struct Best {
        qlen: u64,
        qcovs: f64,
        pident: f64,
        evalue: f64,
        bitscore: f64,
    }
    let mut hits: HashMap<String, Best> = HashMap::new();
    let mut text = String::new();
    std::fs::File::open(&out)?.read_to_string(&mut text)?;
    for line in text.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 8 {
            continue;
        }
        let gene_id = f[0].to_string();
        let pident: f64 = f[2].parse().unwrap_or(0.0);
        let qlen: u64 = f[4].parse().unwrap_or(0);
        let qcovs: f64 = f[5].parse().unwrap_or(0.0);
        let evalue: f64 = f[6].parse().unwrap_or(f64::INFINITY);
        let bitscore: f64 = f[7].parse().unwrap_or(0.0);
        match hits.get_mut(&gene_id) {
            Some(b) if b.bitscore >= bitscore => {}
            _ => {
                hits.insert(
                    gene_id,
                    Best {
                        qlen,
                        qcovs,
                        pident,
                        evalue,
                        bitscore,
                    },
                );
            }
        }
    }

    let mut rows = Vec::new();
    let mut ids: Vec<String> = hits.keys().cloned().collect();
    ids.sort();
    for id in ids {
        let h = &hits[&id];
        // R: qcovs >= blast_cov AND pident >= blast_pid -> Present,
        // any other hit -> Fragment/low identity (PARTIAL here)
        let call = if h.qcovs >= blast_cov && h.pident >= blast_pid {
            Call::Present
        } else {
            Call::Partial
        };
        rows.push(PanelRow {
            gene_id: id,
            qlen: h.qlen,
            cov_pct: h.qcovs,
            identity: h.pident,
            best_evalue: format_evalue(h.evalue),
            call,
        });
    }
    // Panel genes without any hit: read the panel fasta for the id list
    // (and their real lengths).
    let panel_recs = crate::fasta::parse_fasta(panel_fasta)?;
    for rec in panel_recs {
        if !hits.contains_key(&rec.id) {
            rows.push(PanelRow {
                gene_id: rec.id,
                qlen: rec.seq.len() as u64,
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
