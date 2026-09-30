//! Where a panel gene sits in a query genome, and how much of the element
//! carrying it the other genomes of a run hold.
//!
//! The first question is answered from what a run already has: the
//! query contigs, how much of each aligns to the reference (the nucmer
//! delta), and the gained regions with their predicted, named genes. A
//! gene on a contig that barely aligns to the reference chromosome, next
//! to a transposase or a plasmid replication protein, is most likely on
//! a plasmid or a mobile element; the verdict says so as an indication
//! and shows the numbers behind it.
//!
//! The second question takes a whole element (the contig carrying the
//! gene in one query, or a complete plasmid from NCBI) and BLASTs it
//! against every query, so "the gene is missing" can be told apart from
//! "the whole plasmid is missing". Absence from an assembly is still not
//! absence from the bacterium: a small plasmid can be lost in the lab or
//! in library preparation. The interface says so; the engine only counts.

use crate::delta::{merge_intervals, DeltaFile};
use crate::fasta::{self, FastaRecord};
use crate::tools::ToolPaths;
use crate::{friendly, EngineError, Result};
use std::collections::HashMap;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;
use straincompass_types::{ContextGene, ContigStat, ElementHit, GainedAnchor, GainedRow};

/// Genes this far either side of the hit are listed as its context.
pub const CONTEXT_WINDOW: u64 = 10_000;

/// A contig with less than this share aligned to the reference does not
/// look like the reference chromosome.
const PLASMID_MAX_ALIGNED: f64 = 0.2;

/// Listeria plasmids run from ~2 kb to ~100 kb; a contig much larger
/// than any plasmid that still does not align is more likely a divergent
/// piece of chromosome.
const PLASMID_MAX_LEN: u64 = 500_000;

/// A panel hit's place in the query: (contig, start, end, strand).
pub fn parse_locus(locus: &str) -> Option<(String, u64, u64, i8)> {
    let (contig, rest) = locus.rsplit_once(':')?;
    let (range, strand) = match rest.strip_suffix(")") {
        Some(r) => {
            let (range, s) = r.rsplit_once('(')?;
            (range, if s == "-" { -1 } else { 1 })
        }
        None => (rest, 1),
    };
    let (a, b) = range.split_once('-')?;
    Some((contig.to_string(), a.parse().ok()?, b.parse().ok()?, strand))
}

/// Length and reference-aligned bases of every query contig, in fasta order.
pub fn contig_stats(qry_records: &[FastaRecord], delta: &DeltaFile) -> Vec<ContigStat> {
    let aligned = delta.aligned_qry_intervals();
    qry_records
        .iter()
        .map(|r| ContigStat {
            seqid: r.id.clone(),
            length: r.seq.len() as u64,
            aligned_bp: aligned
                .get(&r.id)
                .map(|ivs| ivs.iter().map(|(s, e)| e - s + 1).sum())
                .unwrap_or(0),
        })
        .collect()
}

/// Whether a gene name points at a mobile element.
pub fn is_mobile_name(label: &str) -> bool {
    let l = label.to_lowercase();
    const WORDS: [&str; 12] = [
        "transposase",
        "resolvase",
        "integrase",
        "recombinase",
        "insertion element",
        "insertion sequence",
        "is family",
        "plasmid",
        "replication protein",
        "replication initiat",
        "mobilization",
        "relaxase",
    ];
    if WORDS.iter().any(|w| l.contains(w)) {
        return true;
    }
    // gene symbols: tnpA, tnpR, repA, mobA ...
    let sym = label.trim();
    ["tnp", "rep", "mob"].iter().any(|p| {
        sym.len() >= 4
            && sym.len() <= 6
            && sym.to_lowercase().starts_with(p)
            && sym[3..].chars().all(|c| c.is_ascii_uppercase() || c.is_ascii_digit())
    })
}

/// The predicted genes of the gained regions within `window` of the hit.
pub fn context_genes(
    rows: &[GainedRow],
    contig: &str,
    hit: (u64, u64),
    window: u64,
) -> Vec<ContextGene> {
    let lo = hit.0.saturating_sub(window);
    let hi = hit.1 + window;
    let mut out = Vec::new();
    for r in rows.iter().filter(|r| r.qry_seqid == contig) {
        for o in &r.orfs {
            if o.end < lo || o.start > hi {
                continue;
            }
            let (label, source) = match (&o.best, &o.ncbi) {
                (Some(b), _) if !b.label.is_empty() => (b.label.clone(), "reference"),
                (Some(b), None) => (b.locus_tag.clone(), "reference"),
                (_, Some(n)) => (n.label.clone(), "ncbi"),
                _ => (String::new(), ""),
            };
            let distance = if o.end < hit.0 {
                hit.0 - o.end
            } else {
                // 0 when the gene overlaps the hit
                o.start.saturating_sub(hit.1)
            };
            out.push(ContextGene {
                start: o.start,
                end: o.end,
                strand: o.strand,
                mobile: is_mobile_name(&label),
                label,
                source: source.into(),
                distance,
                is_hit: distance == 0,
                partial: o.partial,
            });
        }
    }
    out.sort_by_key(|g| g.start);
    out
}

/// The verdict for one hit: a short key and the sentence behind it.
pub fn verdict(
    contig: &ContigStat,
    region: Option<&GainedRow>,
    mobile_markers: &[String],
) -> (String, String) {
    let frac = if contig.length > 0 {
        contig.aligned_bp as f64 / contig.length as f64
    } else {
        0.0
    };
    let pct = 100.0 * frac;
    let size = human_bp(contig.length);
    let (key, mut text) = if frac < PLASMID_MAX_ALIGNED && contig.length <= PLASMID_MAX_LEN {
        (
            "plasmid",
            format!(
                "Likely extrachromosomal (plasmid): the gene is on {} ({size}), of which only {pct:.0} % aligns to the reference chromosome.",
                contig.seqid
            ),
        )
    } else if frac < PLASMID_MAX_ALIGNED {
        (
            "unplaced",
            format!(
                "The gene is on {} ({size}), a contig that barely aligns to the reference ({pct:.0} %); it may be a divergent part of the chromosome.",
                contig.seqid
            ),
        )
    } else {
        match region {
            Some(r) if r.anchor == GainedAnchor::Between => (
                "chromosome_insertion",
                format!(
                    "Inserted in the chromosome: the gene lies in a {} stretch with no counterpart in the reference, between {} and {} ({} is {pct:.0} % aligned to the reference).",
                    human_bp(r.length),
                    or_dash(&r.left_gene),
                    or_dash(&r.right_gene),
                    contig.seqid
                ),
            ),
            Some(r) => (
                "chromosome_insertion",
                format!(
                    "Next to chromosomal sequence: the gene lies in a {} stretch with no counterpart in the reference, at the edge of an aligned part of {} ({pct:.0} % aligned); where it is inserted cannot be told from this assembly.",
                    human_bp(r.length),
                    contig.seqid
                ),
            ),
            None => (
                "chromosome_shared",
                format!(
                    "In the chromosome, in a part shared with the reference ({} is {pct:.0} % aligned to the reference).",
                    contig.seqid
                ),
            ),
        }
    };
    if !mobile_markers.is_empty() {
        text.push_str(&format!(
            " Mobile-element genes nearby: {}.",
            mobile_markers.join(", ")
        ));
    }
    (key.into(), text)
}

fn or_dash(s: &str) -> &str {
    if s.is_empty() {
        "-"
    } else {
        s
    }
}

pub fn human_bp(n: u64) -> String {
    if n >= 1_000_000 {
        format!("{:.2} Mb", n as f64 / 1e6)
    } else if n >= 10_000 {
        format!("{:.0} kb", n as f64 / 1e3)
    } else if n >= 1_000 {
        format!("{:.1} kb", n as f64 / 1e3)
    } else {
        format!("{n} bp")
    }
}

/// One genome to search the element in.
pub struct ElementTarget {
    pub query_id: i64,
    pub fasta: PathBuf,
    /// An existing nucleotide BLAST database of `fasta` (the panel
    /// recheck's), reused when its files are still there.
    pub db: Option<PathBuf>,
}

/// Matched stretches of the element in one target, before names are added.
#[derive(Debug, Default, PartialEq)]
pub struct Coverage {
    pub covered_bp: u64,
    pub pieces: usize,
    pub largest_piece: u64,
    pub identity: f64,
    pub contigs: Vec<String>,
}

/// Merge blast hits (element start, element end, identity, aligned
/// length, target contig) into coverage of the element.
pub fn coverage_from_hits(hits: &[(u64, u64, f64, u64, String)]) -> Coverage {
    let mut ivs: Vec<(u64, u64)> = hits.iter().map(|h| (h.0.min(h.1), h.0.max(h.1))).collect();
    merge_intervals(&mut ivs);
    let covered_bp = ivs.iter().map(|(s, e)| e - s + 1).sum();
    let largest_piece = ivs.iter().map(|(s, e)| e - s + 1).max().unwrap_or(0);
    let aligned: u64 = hits.iter().map(|h| h.3).sum();
    let identity = if aligned > 0 {
        hits.iter().map(|h| h.2 * h.3 as f64).sum::<f64>() / aligned as f64
    } else {
        0.0
    };
    let mut by_contig: HashMap<&str, u64> = HashMap::new();
    for h in hits {
        *by_contig.entry(h.4.as_str()).or_default() += h.3;
    }
    let mut contigs: Vec<(&str, u64)> = by_contig.into_iter().collect();
    contigs.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(b.0)));
    Coverage {
        covered_bp,
        pieces: ivs.len(),
        largest_piece,
        identity,
        contigs: contigs.into_iter().map(|(c, _)| c.to_string()).collect(),
    }
}

fn db_ready(db: &Path) -> bool {
    ["nsq", "nal"]
        .iter()
        .any(|ext| PathBuf::from(format!("{}.{ext}", db.display())).is_file())
}

/// BLAST one element against every target genome.
pub fn element_presence(
    tools: &ToolPaths,
    element_id: &str,
    element_seq: &[u8],
    targets: &[ElementTarget],
    work: &Path,
) -> Result<Vec<(i64, Coverage)>> {
    if element_seq.len() < 50 {
        return Err(friendly("The element is too short to compare."));
    }
    std::fs::create_dir_all(work)?;
    let element_fa = work.join("element.fa");
    {
        let mut w = std::io::BufWriter::new(std::fs::File::create(&element_fa)?);
        writeln!(w, ">{}", fasta::sanitize_id(element_id))?;
        for chunk in element_seq.chunks(60) {
            w.write_all(chunk)?;
            writeln!(w)?;
        }
    }
    let mut out = Vec::with_capacity(targets.len());
    for t in targets {
        let db = match &t.db {
            Some(db) if db_ready(db) => db.clone(),
            _ => {
                let db = work.join(format!("db_{}", t.query_id));
                let o = Command::new(&tools.makeblastdb)
                    .arg("-in")
                    .arg(&t.fasta)
                    .args(["-dbtype", "nucl", "-out"])
                    .arg(&db)
                    .output()
                    .map_err(|e| EngineError::ToolMissing(format!("makeblastdb: {e}")))?;
                if !o.status.success() {
                    return Err(friendly(format!(
                        "A genome could not be prepared for the element search. {}",
                        String::from_utf8_lossy(&o.stderr).trim()
                    )));
                }
                db
            }
        };
        let hits_path = work.join(format!("hits_{}.tsv", t.query_id));
        let o = Command::new(&tools.blastn)
            .arg("-query")
            .arg(&element_fa)
            .arg("-db")
            .arg(&db)
            .args([
                "-evalue",
                "1e-20",
                "-max_target_seqs",
                "5000",
                "-outfmt",
                "6 qstart qend pident length sseqid",
                "-out",
            ])
            .arg(&hits_path)
            .output()
            .map_err(|e| EngineError::ToolMissing(format!("blastn: {e}")))?;
        if !o.status.success() {
            return Err(friendly(format!(
                "The element search failed. {}",
                String::from_utf8_lossy(&o.stderr).trim()
            )));
        }
        let text = std::fs::read_to_string(&hits_path)?;
        let hits: Vec<(u64, u64, f64, u64, String)> = text
            .lines()
            .filter_map(|l| {
                let f: Vec<&str> = l.split('\t').collect();
                Some((
                    f.first()?.parse().ok()?,
                    f.get(1)?.parse().ok()?,
                    f.get(2)?.parse().ok()?,
                    f.get(3)?.parse().ok()?,
                    f.get(4)?.to_string(),
                ))
            })
            .collect();
        out.push((t.query_id, coverage_from_hits(&hits)));
    }
    Ok(out)
}

/// Turn per-target coverage into report rows.
pub fn element_hit(
    query_id: i64,
    query_name: String,
    call: Option<straincompass_types::Call>,
    cov: Coverage,
    element_len: u64,
    genome_bp: u64,
) -> ElementHit {
    ElementHit {
        query_id,
        query_name,
        call,
        covered_pct: if element_len > 0 {
            100.0 * cov.covered_bp as f64 / element_len as f64
        } else {
            0.0
        },
        covered_bp: cov.covered_bp,
        pieces: cov.pieces,
        largest_piece: cov.largest_piece,
        identity: cov.identity,
        contigs: cov.contigs,
        genome_bp,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use straincompass_types::{GainedOrf, OrfMatch};

    #[test]
    fn parses_panel_loci() {
        assert_eq!(
            parse_locus("LM226_contig1:2027-2413(-)"),
            Some(("LM226_contig1".into(), 2027, 2413, -1))
        );
        assert_eq!(
            parse_locus("NZ_X:1:5-9(+)"),
            Some(("NZ_X:1".into(), 5, 9, 1))
        );
        assert_eq!(parse_locus("LM226_contig1"), None);
        assert_eq!(parse_locus(""), None);
    }

    #[test]
    fn merges_hits_into_element_coverage() {
        let hits = vec![
            (1, 100, 100.0, 100, "c1".to_string()),
            (90, 200, 98.0, 111, "c1".to_string()),
            (501, 600, 90.0, 100, "c2".to_string()),
        ];
        let c = coverage_from_hits(&hits);
        assert_eq!(c.covered_bp, 300);
        assert_eq!(c.pieces, 2);
        assert_eq!(c.largest_piece, 200);
        assert_eq!(c.contigs, vec!["c1".to_string(), "c2".to_string()]);
        assert!((c.identity - (10000.0 + 98.0 * 111.0 + 9000.0) / 311.0).abs() < 1e-9);
        assert_eq!(coverage_from_hits(&[]), Coverage::default());
    }

    #[test]
    fn contig_stats_counts_aligned_query_bases() {
        let delta = crate::delta::parse_delta_str(
            "/r.fa /q.fa\nNUCMER\n>ref q1 1000 500\n1 100 101 200 0 0 0\n0\n1 50 151 250 0 0 0\n0\n",
        )
        .unwrap();
        let recs = vec![
            FastaRecord { id: "q1".into(), desc: String::new(), seq: vec![b'A'; 500] },
            FastaRecord { id: "q2".into(), desc: String::new(), seq: vec![b'A'; 40] },
        ];
        let s = contig_stats(&recs, &delta);
        assert_eq!((s[0].length, s[0].aligned_bp), (500, 150));
        assert_eq!((s[1].length, s[1].aligned_bp), (40, 0));
    }

    #[test]
    fn recognises_mobile_element_names() {
        for yes in ["transposase", "IS30 family transposase", "tnpA", "tnpR", "repA",
                    "Tn3 family resolvase", "plasmid replication protein", "site-specific integrase"] {
            assert!(is_mobile_name(yes), "{yes}");
        }
        for no in ["cadC", "hly", "internalin A", "replicase-associated protein X", "repressor", ""] {
            assert!(!is_mobile_name(no), "{no}");
        }
    }

    fn orf(start: u64, end: u64, label: &str) -> GainedOrf {
        GainedOrf {
            start,
            end,
            strand: 1,
            partial: false,
            confidence: 99.0,
            best: None,
            ncbi: (!label.is_empty()).then(|| OrfMatch {
                label: label.into(),
                ..Default::default()
            }),
        }
    }

    #[test]
    fn verdict_calls_small_unaligned_contigs_plasmids() {
        let region = GainedRow {
            qry_seqid: "c1".into(),
            start: 1,
            end: 4265,
            length: 4265,
            orfs: vec![orf(100, 900, "Tn3 family resolvase"), orf(2027, 2413, "SMR transporter"), orf(20_000, 20_100, "far")],
            ..Default::default()
        };
        let genes = context_genes(std::slice::from_ref(&region), "c1", (2027, 2413), CONTEXT_WINDOW);
        assert_eq!(genes.len(), 2);
        assert!(genes[0].mobile && !genes[0].is_hit && genes[0].distance == 1127);
        assert!(genes[1].is_hit && !genes[1].mobile);
        let contig = ContigStat { seqid: "c1".into(), length: 4265, aligned_bp: 0 };
        let (k, t) = verdict(&contig, Some(&region), &["Tn3 family resolvase".into()]);
        assert_eq!(k, "plasmid");
        assert!(t.contains("0 %") && t.contains("resolvase"), "{t}");

        let chrom = ContigStat { seqid: "c2".into(), length: 900_000, aligned_bp: 880_000 };
        let placed = GainedRow { anchor: GainedAnchor::Between, left_gene: "lmo1".into(), right_gene: "lmo2".into(), length: 3000, ..Default::default() };
        assert_eq!(verdict(&chrom, Some(&placed), &[]).0, "chromosome_insertion");
        assert_eq!(verdict(&chrom, None, &[]).0, "chromosome_shared");
        let big = ContigStat { seqid: "c3".into(), length: 900_000, aligned_bp: 1000 };
        assert_eq!(verdict(&big, None, &[]).0, "unplaced");
    }
}
