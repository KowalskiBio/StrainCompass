//! Strict BLAST recheck of a gene panel against each query genome, and
//! the reference back-check of gained regions.

use crate::gff::Gene;
use crate::tools::ToolPaths;
use crate::{friendly, Result};
use std::collections::HashMap;
use std::io::Read;
use std::path::Path;
use std::process::Command;
use straincompass_types::{
    Call, GainedBlastHit, GainedIdentify, GainedOrf, GainedRow, GainedVerify, IdentifiedOrf,
    OrfMatch, PanelRow,
};

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
/// The input files a gained-region examination reads: the reference
/// with its annotation, and the query fasta the region is cut from.
pub struct RegionInputs<'a> {
    pub ref_fasta: &'a Path,
    pub ref_gff: &'a Path,
    pub qry_fasta: &'a Path,
}

pub fn gained_verify(
    tools: &ToolPaths,
    inputs: &RegionInputs,
    seqid: &str,
    start: u64,
    end: u64,
    work_dir: &Path,
) -> Result<GainedVerify> {
    let records = crate::fasta::parse_fasta(inputs.qry_fasta)?;
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
        .arg(inputs.ref_fasta)
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
    let (mut tx_hits, tx_note) = if hits.is_empty() {
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
    let ref_records = crate::fasta::parse_fasta(inputs.ref_fasta)?;
    let lcs = crate::longest_match::longest_exact_match(&seq, &ref_records);

    // Say what each hit hits: a row of reference coordinates alone
    // answers nothing until it is joined to the annotation.
    let ref_genes = crate::gff::parse_gff(inputs.ref_gff).unwrap_or_default();
    for h in hits
        .iter_mut()
        .chain(weak_hits.iter_mut())
        .chain(tx_hits.iter_mut().flatten())
    {
        h.genes = overlapping_genes(&ref_genes, &h.ref_seqid, h.ref_start, h.ref_end);
    }

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

/// Reference genes overlapping an interval, formatted for display as
/// "locus_tag (symbol)" - symbol preferred, product as fallback,
/// bare locus tag when the annotation says nothing more.
fn overlapping_genes(genes: &[Gene], seqid: &str, start: u64, end: u64) -> Vec<String> {
    let mut out: Vec<String> = Vec::new();
    for g in genes {
        if g.seqid == seqid && g.start <= end && g.end >= start {
            let label = if !g.symbol.is_empty() {
                format!("{} ({})", g.locus_tag, g.symbol)
            } else if !g.product.is_empty() {
                format!("{} ({})", g.locus_tag, truncate_product(&g.product))
            } else {
                g.locus_tag.clone()
            };
            out.push(label);
        }
    }
    out
}

/// A product can be a whole sentence; the label wants a clause.
fn truncate_product(p: &str) -> String {
    const MAX: usize = 60;
    if p.len() <= MAX {
        p.to_string()
    } else {
        format!(
            "{}...",
            &p[..p
                .char_indices()
                .take(MAX)
                .last()
                .map(|(i, _)| i)
                .unwrap_or(MAX)]
        )
    }
}

/// The standard genetic code over the four unambiguous bases, indexed
/// first base x 16 + second x 4 + third (T=0, C=1, A=2, G=3). Table 11,
/// the bacterial code, differs only in which starts count as M, and
/// that difference does not matter to a similarity search. Partial
/// trailing codons are dropped; ambiguous bases translate to X.
fn translate(seq: &[u8]) -> String {
    const TABLE: [char; 64] = [
        'F', 'F', 'L', 'L', 'S', 'S', 'S', 'S', 'Y', 'Y', '*', '*', 'C', 'C', '*', 'W', //
        'L', 'L', 'L', 'L', 'P', 'P', 'P', 'P', 'H', 'H', 'Q', 'Q', 'R', 'R', 'R', 'R', //
        'I', 'I', 'I', 'M', 'T', 'T', 'T', 'T', 'N', 'N', 'K', 'K', 'S', 'S', 'R', 'R', //
        'V', 'V', 'V', 'V', 'A', 'A', 'A', 'A', 'D', 'D', 'E', 'E', 'G', 'G', 'G', 'G',
    ];
    fn base(b: u8) -> Option<usize> {
        match b.to_ascii_uppercase() {
            b'T' => Some(0),
            b'C' => Some(1),
            b'A' => Some(2),
            b'G' => Some(3),
            _ => None,
        }
    }
    seq.chunks(3)
        .filter(|c| c.len() == 3)
        .map(|c| match (base(c[0]), base(c[1]), base(c[2])) {
            (Some(i), Some(j), Some(k)) => TABLE[i * 16 + j * 4 + k],
            _ => 'X',
        })
        .collect()
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
            genes: Vec::new(),
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

/// Name the genes inside one gained region: each predicted ORF searched,
/// as translated DNA, against the reference's own proteins.
///
/// This is the only name an unannotated query genome can be given
/// locally - a query is a bare draft assembly, so its genes arrive with
/// coordinates and nothing else, and the reference annotation is the one
/// naming resource on the server. It names the ORFs the reference does
/// carry homologs of (interrupted copies, repeat-family members); the
/// ORFs it leaves unnamed are exactly the novel ones, for which the
/// response carries the sequences so the client can link out.
///
/// A gene is only offered for matching when its translation is a clean
/// protein: an internal stop means the reference annotation and the
/// coordinates disagree somewhere, and a corrupted db entry would return
/// to haunt every future search.
/// The reference proteome ready for blastx naming: the CDSs translated
/// from the fasta by the GFF, and the blast database built from them.
/// One entry per coding gene; genes whose translation is not a clean
/// protein are skipped - a corrupted db entry would return to haunt
/// every future search.
pub struct ProteinDb {
    /// locus tag -> (protein accession, display label)
    names: HashMap<String, (String, String)>,
    /// The makeblastdb output prefix to search against.
    db: std::path::PathBuf,
}

/// Build [`ProteinDb`] from the reference's annotation. A GFF without
/// any usable CDSs yields an empty database rather than an error: a
/// reference without coding genes cannot name anything, which is an
/// answer, not a failure.
fn build_protein_db(
    tools: &ToolPaths,
    ref_fasta: &Path,
    ref_gff: &Path,
    work_dir: &Path,
) -> Result<ProteinDb> {
    let ref_records = crate::fasta::parse_fasta(ref_fasta)?;
    let by_seqid: HashMap<&str, &crate::fasta::FastaRecord> =
        ref_records.iter().map(|r| (r.id.as_str(), r)).collect();

    let genes = crate::gff::parse_gff(ref_gff)?;
    let mut proteins: Vec<(String, String, String, String)> = Vec::new();
    for g in &genes {
        if g.protein_id.is_empty() {
            // No CDS child joined in: not a coding gene.
            continue;
        }
        let Some(rec) = by_seqid.get(g.seqid.as_str()) else {
            continue;
        };
        let cds = crate::fasta::subseq(rec, g.start, g.end, g.strand < 0);
        let mut prot = translate(&cds);
        // NCBI GFF3 CDS spans include the terminal stop codon: part of
        // the span, not of the protein.
        if prot.ends_with('*') {
            prot.pop();
        }
        if prot.is_empty() || prot.contains('*') {
            continue;
        }
        proteins.push((
            g.locus_tag.clone(),
            g.protein_id.clone(),
            gene_label(g),
            prot,
        ));
    }

    let out = ProteinDb {
        names: proteins
            .iter()
            .map(|(locus, pid, label, _)| (locus.clone(), (pid.clone(), label.clone())))
            .collect(),
        db: work_dir.join("prot_db"),
    };
    if proteins.is_empty() {
        return Ok(out);
    }

    std::fs::create_dir_all(work_dir)?;
    let prot_fa = work_dir.join("ref_prot.fa");
    {
        use std::io::Write;
        let mut w = std::io::BufWriter::new(std::fs::File::create(&prot_fa)?);
        for (locus, protein_id, label, prot) in &proteins {
            writeln!(w, ">{locus} {protein_id} {label}")?;
            for chunk in prot.as_bytes().chunks(60) {
                w.write_all(chunk)?;
                w.write_all(b"\n")?;
            }
        }
        w.flush()?;
    }
    let make_db = Command::new(&tools.makeblastdb)
        .args(["-in"])
        .arg(&prot_fa)
        .args(["-dbtype", "prot", "-out"])
        .arg(&out.db)
        .output()
        .map_err(|e| crate::EngineError::ToolMissing(format!("makeblastdb: {e}")))?;
    if !make_db.status.success() {
        return Err(friendly(format!(
            "The protein database could not be built. {}",
            String::from_utf8_lossy(&make_db.stderr).trim()
        )));
    }
    Ok(out)
}

/// Search one batch of nucleotide sequences (id, plus-strand sequence)
/// against a [`ProteinDb`] with blastx and return each sequence's best
/// match (highest bitscore), as translated DNA: the frame is found by
/// the search, not assumed by the caller.
fn blastx_best(
    tools: &ToolPaths,
    db: &ProteinDb,
    seqs: &[(String, String)],
    work_dir: &Path,
    out_name: &str,
) -> Result<HashMap<String, OrfMatch>> {
    let mut out: HashMap<String, OrfMatch> = HashMap::new();
    if seqs.is_empty() || db.names.is_empty() {
        return Ok(out);
    }
    let query_fa = work_dir.join(format!("{out_name}.fa"));
    {
        use std::io::Write;
        let mut w = std::io::BufWriter::new(std::fs::File::create(&query_fa)?);
        for (id, seq) in seqs {
            writeln!(w, ">{id}")?;
            for chunk in seq.as_bytes().chunks(60) {
                w.write_all(chunk)?;
                w.write_all(b"\n")?;
            }
        }
        w.flush()?;
    }
    let hits_tsv = work_dir.join(format!("{out_name}.tsv"));
    let blast = Command::new(&tools.blastx)
        .args(["-query"])
        .arg(&query_fa)
        .args(["-db"])
        .arg(&db.db)
        .args([
            "-evalue",
            "1e-5",
            "-outfmt",
            "6 qseqid sseqid pident qcovs evalue bitscore",
            "-out",
        ])
        .arg(&hits_tsv)
        .output()
        .map_err(|e| crate::EngineError::ToolMissing(format!("blastx: {e}")))?;
    if !blast.status.success() {
        return Err(friendly(format!(
            "The translated search did not finish correctly. {}",
            String::from_utf8_lossy(&blast.stderr).trim()
        )));
    }
    let mut text = String::new();
    std::fs::File::open(&hits_tsv)?.read_to_string(&mut text)?;
    for line in text.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 6 {
            continue;
        }
        let better = out
            .get(f[0])
            .map(|m: &OrfMatch| f[5].parse::<f64>().unwrap_or(0.0) > m.bitscore)
            .unwrap_or(true);
        if !better {
            continue;
        }
        let Some((protein_id, label)) = db.names.get(f[1]) else {
            continue;
        };
        out.insert(
            f[0].to_string(),
            OrfMatch {
                locus_tag: f[1].to_string(),
                protein_id: protein_id.clone(),
                label: label.clone(),
                identity: f[2].parse().unwrap_or(0.0),
                coverage: f[3].parse().unwrap_or(0.0),
                evalue: f[4].parse().unwrap_or(f64::INFINITY),
                bitscore: f[5].parse().unwrap_or(0.0),
            },
        );
    }
    Ok(out)
}

/// Name the genes inside one gained region: each predicted ORF searched,
/// as translated DNA, against the reference's own proteins.
///
/// This is the only name an unannotated query genome can be given
/// locally - a query is a bare draft assembly, so its genes arrive with
/// coordinates and nothing else, and the reference annotation is the one
/// naming resource on the server. It names the ORFs the reference does
/// carry homologs of (interrupted copies, repeat-family members); the
/// ORFs it leaves unnamed are exactly the novel ones, for which the
/// response carries the sequences so the client can link out.
pub fn gained_identify(
    tools: &ToolPaths,
    inputs: &RegionInputs,
    seqid: &str,
    start: u64,
    end: u64,
    orfs: &[GainedOrf],
    work_dir: &Path,
) -> Result<GainedIdentify> {
    let qry_records = crate::fasta::parse_fasta(inputs.qry_fasta)?;
    let qry_rec = qry_records.iter().find(|r| r.id == seqid).ok_or_else(|| {
        friendly(format!(
            "The query contig {seqid} is not in this query's fasta file."
        ))
    })?;
    let region_seq = crate::fasta::subseq(qry_rec, start, end, false);

    let mut out = GainedIdentify {
        qry_seqid: seqid.to_string(),
        start,
        end,
        orfs: orfs
            .iter()
            .map(|o| IdentifiedOrf {
                start: o.start,
                end: o.end,
                strand: o.strand,
                best: None,
                seq: String::from_utf8_lossy(&crate::fasta::subseq(
                    qry_rec,
                    o.start,
                    o.end,
                    o.strand < 0,
                ))
                .into_owned(),
            })
            .collect(),
        region_seq: String::from_utf8_lossy(&region_seq).into_owned(),
    };

    let db = build_protein_db(tools, inputs.ref_fasta, inputs.ref_gff, work_dir)?;
    let seqs: Vec<(String, String)> = out
        .orfs
        .iter()
        .enumerate()
        .map(|(i, o)| (format!("o{i}"), o.seq.clone()))
        .collect();
    let best = blastx_best(tools, &db, &seqs, work_dir, "orf_hits")?;
    for (i, orf) in out.orfs.iter_mut().enumerate() {
        if let Some(m) = best.get(&format!("o{i}")) {
            orf.best = Some(m.clone());
        }
    }
    Ok(out)
}

/// Name the predicted genes of every gained region of one query in one
/// translated search, filling each ORF's `best` and each row's
/// `gene_names` in place. Called by the run's pipeline right after
/// prodigal, so the names are in the table from the start instead of
/// behind a per-region button.
///
/// Like gene prediction this is decoration on a result that is already
/// complete without it, so the caller degrades to unnamed genes rather
/// than failing the comparison; errors are reported by the on-demand
/// identification endpoint, which exists precisely because a name can
/// also be requested later.
pub fn name_gained_orfs(
    tools: &ToolPaths,
    ref_fasta: &Path,
    ref_gff: &Path,
    rows: &mut [GainedRow],
    by_id: &HashMap<&str, &crate::fasta::FastaRecord>,
    work_dir: &Path,
) -> Result<()> {
    let db = build_protein_db(tools, ref_fasta, ref_gff, work_dir)?;
    // One id per ORF across all regions: r{region}o{orf}, so a single
    // blastx run names the whole query and the results map straight
    // back onto the rows they came from.
    let mut seqs: Vec<(String, String)> = Vec::new();
    for (ri, row) in rows.iter().enumerate() {
        let Some(rec) = by_id.get(row.qry_seqid.as_str()) else {
            continue;
        };
        for (oi, o) in row.orfs.iter().enumerate() {
            seqs.push((
                format!("r{ri}o{oi}"),
                String::from_utf8_lossy(&crate::fasta::subseq(rec, o.start, o.end, o.strand < 0))
                    .into_owned(),
            ));
        }
    }
    let best = blastx_best(tools, &db, &seqs, work_dir, "gained_names")?;
    for (ri, row) in rows.iter_mut().enumerate() {
        let mut names: Vec<String> = Vec::new();
        for (oi, o) in row.orfs.iter_mut().enumerate() {
            if let Some(m) = best.get(&format!("r{ri}o{oi}")) {
                names.push(if m.label.is_empty() {
                    m.locus_tag.clone()
                } else {
                    m.label.clone()
                });
                o.best = Some(m.clone());
            }
        }
        row.gene_names = names;
        // Marks "the pass ran", so an empty gene_names further down the
        // line reads as genuinely unmatched rather than unknown.
        row.named = true;
    }
    Ok(())
}

/// The display label of a gene: its symbol when annotated, else its
/// product, else nothing (the locus tag carries it alone).
fn gene_label(g: &Gene) -> String {
    if !g.symbol.is_empty() {
        g.symbol.clone()
    } else {
        g.product.clone()
    }
}
