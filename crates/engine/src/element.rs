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
use straincompass_types::{
    Call, ContextGene, ContigStat, ElementHit, GainedAnchor, GainedRow, PanelRow, WgaGene,
    SHORT_MATCH_MAX_COVERAGE, VARIANT_MIN_PROTEIN_IDENTITY,
};

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
            && sym[3..]
                .chars()
                .all(|c| c.is_ascii_uppercase() || c.is_ascii_digit())
    })
}

/// Where a panel gene's sequence came from, read from its record in the
/// panel FASTA as the panel builder writes it: "AMRFinderPlus emrC_Lis:
/// product [origin]", "VFDB ...", "NCBI <accession>", "reference
/// <locus>". A record without such a note is matched against the
/// reference's genes by name, else it is the user's own sequence.
pub fn panel_gene_source(id: &str, desc: &str, ref_genes: &[WgaGene]) -> String {
    let mut words = desc.split_whitespace();
    let first = words.next().unwrap_or("");
    let entry = words.next().unwrap_or("").trim_end_matches(':');
    let with = |db: &str| {
        if entry.is_empty() {
            db.to_string()
        } else {
            format!("{db} ({entry})")
        }
    };
    match first {
        "AMRFinderPlus" | "VFDB" => with(first),
        "NCBI" => with("NCBI Nucleotide"),
        "reference" => with("reference genome"),
        _ => match ref_genes.iter().find(|g| {
            g.locus_tag == id || (!g.symbol.is_empty() && g.symbol.eq_ignore_ascii_case(id))
        }) {
            Some(g) => format!("reference genome ({})", g.locus_tag),
            None => "your panel FASTA".into(),
        },
    }
}

/// Where to show a panel gene in a query: its DNA hit, else the place of
/// its protein-level relative; empty when neither exists. An absent gene
/// with a relative is shown at the relative, the whole gene, rather than
/// at the short stretch the DNA search may have matched inside it.
pub fn panel_locus(row: &PanelRow) -> &str {
    match (&row.protein, row.call) {
        (Some(p), Call::Absent) => &p.locus,
        _ if !row.qry_locus.is_empty() => &row.qry_locus,
        (Some(p), _) => &p.locus,
        _ => "",
    }
}

fn loci_overlap(a: &str, b: &str) -> bool {
    match (parse_locus(a), parse_locus(b)) {
        (Some((ca, sa, ea, _)), Some((cb, sb, eb, _))) => ca == cb && sa <= eb && sb <= ea,
        _ => false,
    }
}

/// Settle each panel call against what the protein search found.
///
/// A protein match is another variant of the gene only when nothing else
/// explains it. It is not when it sits on another panel gene found in the
/// strain (qacH's match is the strain's emrC, a ~70 % relative), nor when
/// it lies in DNA the strain shares with the reference, where it is a
/// chromosomal relative the reference has too (gadD1 against the core
/// glutamate decarboxylase). `gained` holds the strain's stretches with
/// no reference counterpart; None skips that test.
///
/// Then a short partial match (under `SHORT_MATCH_MAX_COVERAGE` of the
/// gene) whose protein match is no variant is a conserved stretch of
/// another gene of the family, and the gene is absent.
pub fn settle_panel_calls(rows: &mut [PanelRow], gained: Option<&[GainedRow]>) {
    let found: Vec<(String, String)> = rows
        .iter()
        .filter(|r| r.call != Call::Absent && !r.qry_locus.is_empty())
        .map(|r| (r.gene_id.clone(), r.qry_locus.clone()))
        .collect();
    for r in rows.iter_mut() {
        let Some(p) = r.protein.as_mut() else {
            continue;
        };
        let other = found
            .iter()
            .find(|(g, l)| *g != r.gene_id && loci_overlap(l, &p.locus));
        if let Some((g, _)) = other {
            p.explained_by = format!("{g}, another panel gene found in this strain");
        } else if let (Some(gained), Some((contig, s, e, _))) = (gained, parse_locus(&p.locus)) {
            let in_gained = gained
                .iter()
                .any(|x| x.qry_seqid == contig && x.start <= e && s <= x.end);
            if !in_gained {
                p.explained_by =
                    "a gene in DNA this strain shares with the reference genome".into();
            }
        }
        r.variant_warning = p.explained_by.is_empty() && p.suggests_variant();
        if r.call == Call::Partial
            && r.cov_pct < SHORT_MATCH_MAX_COVERAGE
            && !r.variant_warning
            && p.identity < VARIANT_MIN_PROTEIN_IDENTITY
        {
            r.call = Call::Absent;
        }
    }
}

/// What a panel hit is when it is not a full match, in words: a partial
/// DNA match, or only a protein-level relative, and whether that relative
/// is close enough to be another variant of the gene. Empty for a gene
/// found in full, or not found at all by a run that checked proteins.
pub fn panel_match_note(row: &PanelRow) -> String {
    let g = &row.gene_id;
    let mut s = match (row.call, &row.protein) {
        (Call::Present, _) => return String::new(),
        (Call::Partial, _) => format!(
            "Only part of {g} matches at DNA level: {:.0} % identity over {:.0} % of the gene.",
            row.identity, row.cov_pct
        ),
        (Call::Absent, Some(p)) if row.cov_pct > 0.0 => format!(
            "{g} is not here. Only a short stretch, {:.0} % of the gene at {:.0} % identity, matches at DNA level, inside a related gene at {}.",
            row.cov_pct, row.identity, p.locus
        ),
        (Call::Absent, Some(p)) => format!(
            "No DNA match for {g}. A protein-level search found a related gene at {}.",
            p.locus
        ),
        (Call::Absent, None) if row.n_variants == 0 => {
            return format!(
                "This run predates the protein-level check; run the comparison again to look for other variants of {g}."
            )
        }
        (Call::Absent, None) => return String::new(),
    };
    if let Some(p) = &row.protein {
        s.push_str(&format!(
            " At protein level it is {:.0} % identical over {:.0} % of the protein.",
            p.identity, p.coverage
        ));
        if !p.explained_by.is_empty() {
            s.push_str(&format!(
                " That match is {}, not another variant of {g}.",
                p.explained_by
            ));
        } else if row.variant_warning {
            s.push_str(&format!(
                " That is close enough to be another variant of {g}, one the panel does not hold: check it before calling {g} absent."
            ));
        } else {
            s.push_str(&format!(
                " That is too distant to be {g} itself; most likely a related gene from the same family."
            ));
        }
    }
    s
}

/// Share of both the predicted gene and the hit that must overlap for the
/// gene to count as the panel gene itself, percent.
const SAME_LOCUS_MIN_OVERLAP: f64 = 80.0;

/// The two intervals overlap over at least `SAME_LOCUS_MIN_OVERLAP` % of
/// each, so they are the same gene seen twice.
fn same_locus(a: (u64, u64), b: (u64, u64)) -> bool {
    let shared = a.1.min(b.1) + 1;
    let shared = shared.saturating_sub(a.0.max(b.0)) as f64;
    let len = |x: (u64, u64)| (x.1 + 1).saturating_sub(x.0).max(1) as f64;
    shared * 100.0 >= SAME_LOCUS_MIN_OVERLAP * len(a)
        && shared * 100.0 >= SAME_LOCUS_MIN_OVERLAP * len(b)
}

/// The predicted genes of the gained regions within `window` of the hit.
///
/// The predicted gene that is the hit is named after the panel gene
/// (`gene_id`): the panel search has already identified it, while its
/// reference-proteome name is only the closest relative the reference
/// has (an emrC hit came out as the reference's "sugE2"). That name is
/// kept as `match_label`, with its identity, so the reader can judge it.
pub fn context_genes(
    rows: &[GainedRow],
    contig: &str,
    hit: (u64, u64),
    window: u64,
    gene_id: &str,
) -> Vec<ContextGene> {
    let lo = hit.0.saturating_sub(window);
    let hi = hit.1 + window;
    let mut out = Vec::new();
    for r in rows.iter().filter(|r| r.qry_seqid == contig) {
        for o in &r.orfs {
            if o.end < lo || o.start > hi {
                continue;
            }
            let (m, source) = match (&o.best, &o.ncbi) {
                (Some(b), _) if !b.label.is_empty() || o.ncbi.is_none() => (Some(b), "reference"),
                (_, Some(n)) => (Some(n), "ncbi"),
                _ => (None, ""),
            };
            let name = m
                .map(|m| {
                    if m.label.is_empty() {
                        m.locus_tag.clone()
                    } else {
                        m.label.clone()
                    }
                })
                .unwrap_or_default();
            let distance = if o.end < hit.0 {
                hit.0 - o.end
            } else {
                // 0 when the gene overlaps the hit
                o.start.saturating_sub(hit.1)
            };
            let is_panel = !gene_id.is_empty() && same_locus((o.start, o.end), hit);
            let (label, source, match_label) = if is_panel {
                (gene_id.to_string(), "panel", name)
            } else {
                (name, source, String::new())
            };
            out.push(ContextGene {
                start: o.start,
                end: o.end,
                strand: o.strand,
                mobile: is_mobile_name(&label),
                label,
                source: source.into(),
                match_label,
                match_identity: m.map(|m| m.identity),
                match_coverage: m.map(|m| m.coverage).filter(|c| *c > 0.0),
                locus_tag: String::new(),
                distance,
                is_hit: distance == 0,
                partial: o.partial,
            });
        }
    }
    out.sort_by_key(|g| g.start);
    out
}

/// The reference's own genes in the parts of the window that align to
/// the reference, placed on the query contig through the alignment
/// blocks. Predicted genes exist only in gained regions, so without these
/// a gene in a shared stretch (plcA in LIPI-1) has no neighbours at all.
///
/// Positions are mapped by offset within each block, which ignores the
/// block's few indels: good to a handful of bases, plenty for a list.
pub fn annotation_context_genes(
    alignments: &[crate::delta::Alignment],
    genes: &[WgaGene],
    contig: &str,
    hit: (u64, u64),
    window: u64,
) -> Vec<ContextGene> {
    let lo = hit.0.saturating_sub(window);
    let hi = hit.1 + window;
    let mut by_seqid: HashMap<&str, Vec<&WgaGene>> = HashMap::new();
    for g in genes {
        by_seqid.entry(g.seqid.as_str()).or_default().push(g);
    }
    let mut out: Vec<ContextGene> = Vec::new();
    for a in alignments.iter().filter(|a| a.qry_seqid == contig) {
        let (qs, qe) = (a.qry_lo.max(lo), a.qry_hi.min(hi));
        if qs > qe {
            continue;
        }
        // the clipped query stretch on the reference
        let to_ref = |q: u64| {
            if a.qry_rev {
                a.ref_end.saturating_sub(q - a.qry_lo)
            } else {
                a.ref_start + (q - a.qry_lo)
            }
        };
        let (rs, re) = {
            let (x, y) = (to_ref(qs), to_ref(qe));
            (x.min(y).max(a.ref_start), x.max(y).min(a.ref_end))
        };
        let to_qry = |r: u64| {
            if a.qry_rev {
                a.qry_lo + (a.ref_end - r)
            } else {
                a.qry_lo + (r - a.ref_start)
            }
        };
        let Some(cands) = by_seqid.get(a.ref_seqid.as_str()) else {
            continue;
        };
        for g in cands.iter().filter(|g| g.end >= rs && g.start <= re) {
            // the part of the gene this block covers, back on the query
            let (gs, ge) = (g.start.max(a.ref_start), g.end.min(a.ref_end));
            let (x, y) = (to_qry(gs), to_qry(ge));
            let (start, end) = (x.min(y), x.max(y));
            let distance = if end < hit.0 {
                hit.0 - end
            } else {
                start.saturating_sub(hit.1)
            };
            let product = gff_unescape(&g.product);
            let label = if !g.symbol.is_empty() {
                g.symbol.clone()
            } else {
                product.clone()
            };
            let gene = ContextGene {
                start,
                end,
                strand: if a.qry_rev { -g.strand } else { g.strand },
                mobile: is_mobile_name(&label) || is_mobile_name(&product),
                label,
                source: "annotation".into(),
                locus_tag: g.locus_tag.clone(),
                distance,
                is_hit: distance == 0,
                partial: g.start < a.ref_start || g.end > a.ref_end,
                ..Default::default()
            };
            // a gene split over two blocks, or seen through a repeat
            // block too: keep the copy closest to the hit
            match out.iter_mut().find(|o| o.locus_tag == gene.locus_tag) {
                Some(o) if o.distance > gene.distance => *o = gene,
                Some(_) => {}
                None => out.push(gene),
            }
        }
    }
    out.sort_by_key(|g| g.start);
    out
}

/// Undo GFF3 attribute escaping ("ABC transporter%2C ATP-binding").
fn gff_unescape(v: &str) -> String {
    let b = v.as_bytes();
    let mut out = Vec::with_capacity(b.len());
    let mut i = 0;
    while i < b.len() {
        if b[i] == b'%'
            && i + 2 < b.len()
            && b[i + 1].is_ascii_hexdigit()
            && b[i + 2].is_ascii_hexdigit()
        {
            let hex = |c: u8| (c as char).to_digit(16).unwrap_or(0) as u8;
            out.push(hex(b[i + 1]) * 16 + hex(b[i + 2]));
            i += 3;
            continue;
        }
        out.push(b[i]);
        i += 1;
    }
    String::from_utf8(out).unwrap_or_else(|_| v.to_string())
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
            FastaRecord {
                id: "q1".into(),
                desc: String::new(),
                seq: vec![b'A'; 500],
            },
            FastaRecord {
                id: "q2".into(),
                desc: String::new(),
                seq: vec![b'A'; 40],
            },
        ];
        let s = contig_stats(&recs, &delta);
        assert_eq!((s[0].length, s[0].aligned_bp), (500, 150));
        assert_eq!((s[1].length, s[1].aligned_bp), (40, 0));
    }

    #[test]
    fn places_reference_genes_through_the_alignment() {
        // query 1001-3000 = ref 501-2500 forward; query 5001-6000 = ref
        // 9001-10000 reversed
        let delta = crate::delta::parse_delta_str(
            "/r.fa /q.fa\nNUCMER\n>ref q1 20000 9000\n501 2500 1001 3000 0 0 0\n0\n9001 10000 6000 5001 0 0 0\n0\n",
        )
        .unwrap();
        let g = |tag: &str, sym: &str, s: u64, e: u64, strand: i8| WgaGene {
            locus_tag: tag.into(),
            symbol: sym.into(),
            biotype: "protein_coding".into(),
            seqid: "ref".into(),
            start: s,
            end: e,
            strand,
            product: "transposase".into(),
        };
        let genes = vec![
            g("t1", "plcA", 1001, 1500, 1),
            g("t2", "", 2400, 2700, 1),
            g("t3", "hly", 9101, 9200, 1),
            g("t4", "far", 15000, 15100, 1),
        ];
        let out = annotation_context_genes(&delta.alignments, &genes, "q1", (1501, 2000), 10_000);
        let tags: Vec<&str> = out.iter().map(|o| o.locus_tag.as_str()).collect();
        assert_eq!(tags, vec!["t1", "t2", "t3"]);
        assert!(out[0].is_hit && !out[0].partial && (out[0].start, out[0].end) == (1501, 2000));
        // clipped at the block end, unnamed so labelled by its product
        assert!(
            out[1].partial && out[1].end == 3000 && out[1].label == "transposase" && out[1].mobile
        );
        // reversed block: ref 9101-9200 -> query 5801-5900 on the other strand
        assert_eq!((out[2].start, out[2].end, out[2].strand), (5801, 5900, -1));
        assert_eq!(out[2].distance, 5801 - 2000);
    }

    #[test]
    fn unescapes_gff_attribute_values() {
        assert_eq!(
            gff_unescape("ABC transporter%2C ATP-binding"),
            "ABC transporter, ATP-binding"
        );
        assert_eq!(gff_unescape("a%3Bb%3Dc%25"), "a;b=c%");
        assert_eq!(gff_unescape("100%"), "100%");
        assert_eq!(gff_unescape("%zz"), "%zz");
        assert_eq!(gff_unescape("%\u{e9}t\u{e9}"), "%\u{e9}t\u{e9}");
    }

    #[test]
    fn recognises_mobile_element_names() {
        for yes in [
            "transposase",
            "IS30 family transposase",
            "tnpA",
            "tnpR",
            "repA",
            "Tn3 family resolvase",
            "plasmid replication protein",
            "site-specific integrase",
        ] {
            assert!(is_mobile_name(yes), "{yes}");
        }
        for no in [
            "cadC",
            "hly",
            "internalin A",
            "replicase-associated protein X",
            "repressor",
            "",
        ] {
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
    fn settles_calls_against_the_protein_search() {
        use straincompass_types::ProteinHit;
        let prot = |identity: f64, locus: &str| {
            Some(ProteinHit {
                variant: String::new(),
                identity,
                coverage: 95.0,
                locus: locus.into(),
                ..Default::default()
            })
        };
        let row = |gene: &str, call: Call, cov: f64, locus: &str| PanelRow {
            gene_id: gene.into(),
            call,
            cov_pct: cov,
            identity: 75.0,
            qry_locus: locus.into(),
            n_variants: 1,
            ..Default::default()
        };
        let mut rows = vec![
            // a 135 bp stretch of another metal pump, 37 % as protein
            PanelRow {
                protein: prot(37.0, "chr:1000-2900(+)"),
                ..row("cadA", Call::Partial, 6.0, "chr:2400-2535(+)")
            },
            // qacH's match is the strain's emrC
            PanelRow {
                protein: prot(71.0, "p1:2087-2413(-)"),
                ..row("qacH", Call::Partial, 87.0, "p1:2087-2413(-)")
            },
            row("emrC", Call::Present, 100.0, "p1:2027-2413(-)"),
            // gadD1 against the core glutamate decarboxylase, shared DNA
            PanelRow {
                protein: prot(70.0, "chr:58314-59516(+)"),
                ..row("lmo0447", Call::Partial, 87.0, "chr:58314-59516(+)")
            },
            // an unexplained 70 % relative on a plasmid: a real candidate
            PanelRow {
                protein: prot(70.0, "p2:300-2400(-)"),
                ..row("cadX", Call::Absent, 0.0, "")
            },
        ];
        let gained = vec![
            GainedRow {
                qry_seqid: "p1".into(),
                start: 1,
                end: 4265,
                ..Default::default()
            },
            GainedRow {
                qry_seqid: "p2".into(),
                start: 1,
                end: 13000,
                ..Default::default()
            },
        ];
        settle_panel_calls(&mut rows, Some(&gained));
        assert_eq!(rows[0].call, Call::Absent);
        assert_eq!(panel_locus(&rows[0]), "chr:1000-2900(+)");
        assert!(panel_match_note(&rows[0]).contains("short stretch"));
        assert!(!rows[1].variant_warning);
        assert!(rows[1]
            .protein
            .as_ref()
            .unwrap()
            .explained_by
            .starts_with("emrC"));
        assert_eq!(rows[1].call, Call::Partial);
        assert!(!rows[3].variant_warning);
        assert!(panel_match_note(&rows[3]).contains("shares with the reference"));
        assert!(rows[4].variant_warning);
        assert!(panel_match_note(&rows[4]).contains("another variant of cadX"));
    }

    #[test]
    fn explains_a_gene_not_found_in_full() {
        use straincompass_types::ProteinHit;
        let tn5422 = ProteinHit {
            variant: "cadA".into(),
            identity: 70.2,
            coverage: 98.0,
            locus: "c46:328-2463(-)".into(),
            ..Default::default()
        };
        let mut row = PanelRow {
            gene_id: "cadA".into(),
            call: Call::Absent,
            n_variants: 1,
            protein: Some(tn5422.clone()),
            variant_warning: tn5422.suggests_variant(),
            ..Default::default()
        };
        assert_eq!(panel_locus(&row), "c46:328-2463(-)");
        let n = panel_match_note(&row);
        assert!(
            n.contains("No DNA match") && n.contains("70 %") && n.contains("another variant"),
            "{n}"
        );

        // a family relative, not the gene
        row.protein = Some(ProteinHit {
            identity: 41.0,
            ..tn5422
        });
        row.variant_warning = false;
        assert!(panel_match_note(&row).contains("same family"));

        // the 12 strains without Tn5422: a 135 bp stretch of another ATPase
        row.call = Call::Partial;
        row.identity = 75.4;
        row.cov_pct = 6.0;
        let n = panel_match_note(&row);
        assert!(
            n.contains("6 % of the gene") && n.contains("same family"),
            "{n}"
        );
        row.call = Call::Absent;

        row.protein = None;
        assert_eq!(panel_match_note(&row), "");
        row.n_variants = 0;
        assert!(panel_match_note(&row).contains("predates"));

        row.call = Call::Present;
        assert_eq!(panel_match_note(&row), "");
    }

    #[test]
    fn names_the_database_a_panel_gene_came_from() {
        let refs = vec![WgaGene {
            locus_tag: "lmo0200".into(),
            symbol: "prfA".into(),
            ..Default::default()
        }];
        let src = |id, desc| panel_gene_source(id, desc, &refs);
        assert_eq!(
            src("emrC", "AMRFinderPlus emrC_Lis: EmrC [EAC4468893.1]"),
            "AMRFinderPlus (emrC_Lis)"
        );
        assert_eq!(src("llsA", "VFDB VFG045 : LLS"), "VFDB (VFG045)");
        assert_eq!(
            src("qacH", "NCBI HF565366.1:1-381"),
            "NCBI Nucleotide (HF565366.1:1-381)"
        );
        assert_eq!(
            src("lmo0200", "reference lmo0200"),
            "reference genome (lmo0200)"
        );
        // panels built before the note was written
        assert_eq!(src("prfA", ""), "reference genome (lmo0200)");
        assert_eq!(src("myGene", "whatever"), "your panel FASTA");
    }

    #[test]
    fn hit_gene_takes_the_panel_name_and_keeps_its_match() {
        let mut hit = orf(2027, 2413, "");
        hit.best = Some(OrfMatch {
            label: "sugE2".into(),
            identity: 41.2,
            coverage: 80.0,
            ..Default::default()
        });
        let mut near = orf(2400, 2900, "");
        near.best = Some(OrfMatch {
            label: "lmo0855".into(),
            identity: 99.0,
            coverage: 100.0,
            ..Default::default()
        });
        let region = GainedRow {
            qry_seqid: "c1".into(),
            orfs: vec![hit, near],
            ..Default::default()
        };
        let genes = context_genes(
            std::slice::from_ref(&region),
            "c1",
            (2030, 2413),
            CONTEXT_WINDOW,
            "emrC",
        );
        let g = &genes[0];
        assert_eq!((g.label.as_str(), g.source.as_str()), ("emrC", "panel"));
        assert_eq!(g.match_label, "sugE2");
        assert_eq!(
            (g.match_identity, g.match_coverage),
            (Some(41.2), Some(80.0))
        );
        // overlapping the hit by a few bases does not make it the hit gene
        let g = &genes[1];
        assert!(g.is_hit);
        assert_eq!(
            (g.label.as_str(), g.source.as_str()),
            ("lmo0855", "reference")
        );
        assert!(g.match_label.is_empty());
    }

    #[test]
    fn verdict_calls_small_unaligned_contigs_plasmids() {
        let region = GainedRow {
            qry_seqid: "c1".into(),
            start: 1,
            end: 4265,
            length: 4265,
            orfs: vec![
                orf(100, 900, "Tn3 family resolvase"),
                orf(2027, 2413, "SMR transporter"),
                orf(20_000, 20_100, "far"),
            ],
            ..Default::default()
        };
        let genes = context_genes(
            std::slice::from_ref(&region),
            "c1",
            (2027, 2413),
            CONTEXT_WINDOW,
            "emrC",
        );
        assert_eq!(genes.len(), 2);
        assert!(genes[0].mobile && !genes[0].is_hit && genes[0].distance == 1127);
        assert!(genes[1].is_hit && !genes[1].mobile);
        let contig = ContigStat {
            seqid: "c1".into(),
            length: 4265,
            aligned_bp: 0,
        };
        let (k, t) = verdict(&contig, Some(&region), &["Tn3 family resolvase".into()]);
        assert_eq!(k, "plasmid");
        assert!(t.contains("0 %") && t.contains("resolvase"), "{t}");

        let chrom = ContigStat {
            seqid: "c2".into(),
            length: 900_000,
            aligned_bp: 880_000,
        };
        let placed = GainedRow {
            anchor: GainedAnchor::Between,
            left_gene: "lmo1".into(),
            right_gene: "lmo2".into(),
            length: 3000,
            ..Default::default()
        };
        assert_eq!(
            verdict(&chrom, Some(&placed), &[]).0,
            "chromosome_insertion"
        );
        assert_eq!(verdict(&chrom, None, &[]).0, "chromosome_shared");
        let big = ContigStat {
            seqid: "c3".into(),
            length: 900_000,
            aligned_bp: 1000,
        };
        assert_eq!(verdict(&big, None, &[]).0, "unplaced");
    }
}
