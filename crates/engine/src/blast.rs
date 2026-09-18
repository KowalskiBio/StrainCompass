//! Strict BLAST recheck of a gene panel against each query genome, and
//! the reference back-check of gained regions.

use crate::tools::ToolPaths;
use crate::{friendly, Result};
use std::collections::HashMap;
use std::io::Read;
use std::path::Path;
use std::process::Command;
use straincompass_types::{Call, GainedBlastHit, GainedVerify, PanelRow};

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

/// The strong tier's E-value ceiling; hits at or below it are what the
/// verdicts are built on.
const STRONG_EVALUE: f64 = 1e-5;
/// The loose ceiling of the whole nucleotide search. Matches between the
/// two ceilings are reported as the weak tier rather than silently
/// dropped.
const WEAK_EVALUE: f64 = 10.0;
/// The translated search runs at most on regions this long: tblastx
/// translates both sides and is easily an order of magnitude slower
/// than blastn, and an unanchored whole contig (a plasmid) would hold a
/// cpu slot for minutes on the shared server.
const TX_MAX_BP: u64 = 50_000;
/// Per-tier hit cap: a region that really is a repeat can produce
/// thousands of HSPs, and the table only ever shows the top.
const MAX_HITS: usize = 50;

/// Back-check one gained region: search its sequence against the whole
/// reference genome.
///
/// The gained table calls a region gained when the whole-genome aligner
/// produced no alignment covering it, and nucmer's default matcher can
/// only seed alignments from matches that are unique on the reference
/// side, so a query copy of a multicopy reference family lands there
/// with no alignment at all. These searches are the closer test of
/// "absent from the reference", and they must therefore be *more*
/// sensitive than the aligner, not equally sensitive: `-task blastn`
/// seeds on 11-mers where the default megablast needs 28, and `-dust no`
/// leaves low-complexity sequence unmasked. For proving absence, every
/// hit counts.
///
/// Three answers come back. The nucleotide search runs once at the
/// loose ceiling and is split into a strong and a weak tier. When the
/// strong tier is empty, a translated search (tblastx) adds the second
/// opinion a divergent coding homolog needs - nucleotide seeding misses
/// relatives below ~70% identity that amino-acid similarity still
/// finds. And independent of both, the longest exact match to the
/// reference is computed cutoff-free (see [`crate::longest_match`]).
///
/// `work_dir` receives the scratch files (the region fasta, the blast
/// database, the hit tables); it must be writable and private to this
/// call. The caller owns its lifetime.
pub fn gained_verify(
    tools: &ToolPaths,
    ref_fasta: &Path,
    qry_fasta: &Path,
    seqid: &str,
    start: u64,
    end: u64,
    work_dir: &Path,
) -> Result<GainedVerify> {
    let records = crate::fasta::parse_fasta(qry_fasta)?;
    let rec = records.iter().find(|r| r.id == seqid).ok_or_else(|| {
        friendly(format!(
            "The query contig {seqid} is not in this query's fasta file."
        ))
    })?;
    let seq = crate::fasta::subseq(rec, start, end, false);
    if seq.is_empty() {
        return Err(friendly(
            "The region is empty or lies outside its query contig.",
        ));
    }

    std::fs::create_dir_all(work_dir)?;
    let region_fasta = work_dir.join("region.fa");
    {
        use std::io::Write;
        let mut w = std::io::BufWriter::new(std::fs::File::create(&region_fasta)?);
        writeln!(w, ">region {seqid}:{start}-{end}")?;
        for chunk in seq.chunks(60) {
            w.write_all(chunk)?;
            w.write_all(b"\n")?;
        }
        w.flush()?;
    }

    let db = work_dir.join("ref_db");
    let make_db = Command::new(&tools.makeblastdb)
        .args(["-in"])
        .arg(ref_fasta)
        .args(["-dbtype", "nucl", "-out"])
        .arg(&db)
        .output()
        .map_err(|e| crate::EngineError::ToolMissing(format!("makeblastdb: {e}")))?;
    if !make_db.status.success() {
        return Err(friendly(format!(
            "The reference back-check could not be prepared. {}",
            String::from_utf8_lossy(&make_db.stderr).trim()
        )));
    }

    // One loose nucleotide search, split into the two tiers afterwards:
    // the strong tier's ceiling is a reporting decision, not a search
    // parameter, so no match inside the loose ceiling is hidden from the
    // scientist who has to believe the verdict.
    let nuc_hits = run_search(
        &tools.blastn,
        &region_fasta,
        &db,
        WEAK_EVALUE,
        work_dir.join("hits.tsv"),
        false,
        "blastn",
    )?;
    let mut hits: Vec<GainedBlastHit> = nuc_hits
        .iter()
        .filter(|h| h.evalue <= STRONG_EVALUE)
        .cloned()
        .collect();
    let mut weak_hits: Vec<GainedBlastHit> = nuc_hits
        .iter()
        .filter(|h| h.evalue > STRONG_EVALUE)
        .cloned()
        .collect();
    top(&mut hits);
    top(&mut weak_hits);

    // The translated search is the second opinion on "found nothing",
    // so it only runs when there is nothing strong to second-guess.
    // tblastx reports its coordinates in nucleotide bases but its
    // identity and length over aligned amino acids.
    let (tx_hits, tx_note) = if hits.is_empty() {
        if seq.len() as u64 <= TX_MAX_BP {
            (
                Some(top_owned(run_search(
                    &tools.tblastx,
                    &region_fasta,
                    &db,
                    STRONG_EVALUE,
                    work_dir.join("tx_hits.tsv"),
                    true,
                    "tblastx",
                )?)),
                None,
            )
        } else {
            (
                None,
                Some(format!(
                    "The translated search was skipped: at {} bp this region is longer than the {} kb cap, and tblastx on it would hold the server for minutes.",
                    seq.len(),
                    TX_MAX_BP / 1000
                )),
            )
        }
    } else {
        (None, None)
    };

    // Cutoff-free and tool-free: the longest exact match anywhere. The
    // reference is parsed here because no caller has it loaded.
    let ref_records = crate::fasta::parse_fasta(ref_fasta)?;
    let lcs = crate::longest_match::longest_exact_match(&seq, &ref_records);

    Ok(GainedVerify {
        qry_seqid: seqid.to_string(),
        start,
        end,
        hits,
        weak_hits,
        tx_hits,
        tx_note,
        longest_exact_bp: lcs.length,
        longest_exact_seqid: lcs.ref_seqid,
        longest_exact_start: lcs.ref_start,
        longest_exact_end: lcs.ref_end,
        longest_exact_qry_start: lcs.qry_start,
        longest_exact_qry_end: lcs.qry_end,
    })
}

/// Run one blast search of the region against the database and parse its
/// tabular output. `translated` picks the tool-appropriate sensitivity
/// flags: blastn seeds on 11-mers with low-complexity unmasked, tblastx
/// (which has no `-task` and filters with SEG on the amino-acid side
/// instead of DUST) gets `-seg no` for the same reason.
fn run_search(
    tool: &Path,
    region_fasta: &Path,
    db: &Path,
    evalue: f64,
    out: std::path::PathBuf,
    translated: bool,
    label: &str,
) -> Result<Vec<GainedBlastHit>> {
    let mut cmd = Command::new(tool);
    cmd.args(["-query"]).arg(region_fasta).args(["-db"]).arg(db);
    if translated {
        cmd.args(["-seg", "no"]);
    } else {
        cmd.args(["-task", "blastn", "-dust", "no"]);
    }
    let blast = cmd
        .args([
            "-evalue",
            &format!("{}", evalue.max(1e-300)),
            "-outfmt",
            "6 sseqid sstart send pident length qstart qend evalue bitscore",
        ])
        .arg("-out")
        .arg(&out)
        .output()
        .map_err(|e| crate::EngineError::ToolMissing(format!("{label}: {e}")))?;
    if !blast.status.success() {
        return Err(friendly(format!(
            "The reference back-check ({label}) did not finish correctly. {}",
            String::from_utf8_lossy(&blast.stderr).trim()
        )));
    }

    let mut hits = Vec::new();
    let mut text = String::new();
    std::fs::File::open(&out)?.read_to_string(&mut text)?;
    for line in text.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 9 {
            continue;
        }
        let (Ok(sstart), Ok(send)) = (f[1].parse::<u64>(), f[2].parse::<u64>()) else {
            continue;
        };
        hits.push(GainedBlastHit {
            ref_seqid: f[0].to_string(),
            // The subject span is reported against the plus strand, so a
            // hit on the reverse strand arrives as sstart > send; the
            // row reports the interval either way, like every other
            // reference coordinate in the app.
            ref_start: sstart.min(send),
            ref_end: sstart.max(send),
            identity: f[3].parse().unwrap_or(0.0),
            length: f[4].parse().unwrap_or(0),
            qry_start: f[5].parse().unwrap_or(1),
            qry_end: f[6].parse().unwrap_or(0),
            evalue: f[7].parse().unwrap_or(f64::INFINITY),
            bitscore: f[8].parse().unwrap_or(0.0),
        });
    }
    Ok(hits)
}

/// Best first, and bounded to `MAX_HITS`.
fn top(hits: &mut Vec<GainedBlastHit>) {
    hits.sort_by(|a, b| b.bitscore.total_cmp(&a.bitscore));
    hits.truncate(MAX_HITS);
}

/// [`top`] for a freshly owned vector, so the call sites stay flat.
fn top_owned(mut hits: Vec<GainedBlastHit>) -> Vec<GainedBlastHit> {
    top(&mut hits);
    hits
}
