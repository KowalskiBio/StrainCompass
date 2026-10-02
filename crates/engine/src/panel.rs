//! Build a gene panel FASTA from a list of gene identifiers.
//!
//! Mirrors how `genes_of_interest.fasta` was produced for the R pipeline:
//! the user supplies identifiers (locus tags like lmo0444, gene symbols
//! like inlA, or "symbol (locus_tag)" pairs), and each gene's sequence is
//! extracted from the reference genome in the gene's own orientation.

/// Marks the extra variants of a panel gene in its FASTA ids: "cadA" is
/// the first sequence for cadA, "cadA__v2" the second. The search counts
/// them as one gene and reports which one matched.
pub const VARIANT_MARK: &str = "__v";

/// The gene a panel record belongs to: "cadA__v2" -> "cadA".
pub fn variant_gene(id: &str) -> &str {
    match id.rsplit_once(VARIANT_MARK) {
        Some((g, n)) if !g.is_empty() && !n.is_empty() && n.bytes().all(|b| b.is_ascii_digit()) => {
            g
        }
        _ => id,
    }
}

/// The FASTA id of the `n`th (1-based) sequence of `gene`.
pub fn variant_id(gene: &str, n: usize) -> String {
    if n <= 1 {
        gene.to_string()
    } else {
        format!("{gene}{VARIANT_MARK}{n}")
    }
}

use crate::fasta::{parse_fasta, revcomp};
use crate::gff::parse_gff;
pub use crate::gff::Gene;
use crate::Result;
use std::collections::HashMap;
use std::path::Path;

pub struct PanelFromIds {
    /// The generated panel FASTA text.
    pub fasta: String,
    /// Identifiers that could not be found in the reference annotation.
    pub missing: Vec<String>,
    /// Identifiers that were resolved, with the panel header used.
    pub found: Vec<String>,
    /// Names the reference gives to more than one gene, in words: which
    /// gene was taken, why, and how to ask for the other.
    pub ambiguous: Vec<String>,
}

/// One reference gene a name could mean, for [`panel_from_ids_with`]'s
/// chooser.
pub struct Candidate<'a> {
    pub gene: &'a Gene,
    /// DNA in the gene's own orientation.
    pub seq: Vec<u8>,
}

/// Picks among the genes of one name: (index, why), or None to take the
/// first in genome order.
pub type Chooser<'c> = dyn Fn(&str, &[Candidate]) -> Option<(usize, String)> + 'c;

/// Extract gene sequences from the reference.
///
/// `ids_text` is the raw content of the user's CSV/TSV/list file; every
/// line holds one gene (extra columns and comments are tolerated, and the
/// parenthesized locus tag in entries like "pva (lmo0446)" wins over the
/// symbol).
pub fn panel_from_ids(ref_fasta: &Path, ref_gff: &Path, ids_text: &str) -> Result<PanelFromIds> {
    panel_from_ids_with(ref_fasta, ref_gff, ids_text, &|_, _| None)
}

/// [`panel_from_ids`] with a say in which gene an ambiguous name means.
/// A name can belong to unrelated genes of one genome: EGD-e calls both
/// the virulence regulator PrfA (lmo0200) and peptide chain release
/// factor 1 (lmo2543) "prfA". `choose` picks one; without an answer the
/// first in genome order is taken. Either way the note says so.
pub fn panel_from_ids_with(
    ref_fasta: &Path,
    ref_gff: &Path,
    ids_text: &str,
    choose: &Chooser,
) -> Result<PanelFromIds> {
    // tolerate a UTF-8 BOM (Excel-style CSVs)
    let ids_text = ids_text.trim_start_matches('\u{feff}');
    let genes = parse_gff(ref_gff)?;
    if genes.is_empty() {
        return Err(crate::friendly(
            "The reference annotation has no genes, so the panel cannot be built from it.",
        ));
    }
    let mut by_locus: HashMap<&str, &Gene> = HashMap::new();
    let mut by_old_locus: HashMap<&str, &Gene> = HashMap::new();
    // every gene of a name, in genome order (a HashMap of one gene per
    // name silently kept the last: "prfA" meant RF1, not PrfA)
    let mut by_symbol: HashMap<String, Vec<&Gene>> = HashMap::new();
    for g in &genes {
        by_locus.insert(g.locus_tag.as_str(), g);
        if !g.old_locus_tag.is_empty() {
            by_old_locus.insert(g.old_locus_tag.as_str(), g);
        }
        if !g.symbol.is_empty() {
            by_symbol
                .entry(g.symbol.to_lowercase())
                .or_default()
                .push(g);
        }
    }
    let mut ambiguous = Vec::new();

    let mut fasta = String::new();
    let mut missing = Vec::new();
    let mut found = Vec::new();
    let mut used: Vec<String> = Vec::new();
    let mut any_entry = false;

    // Lines hold either one gene (one-per-line list/CSV) or several
    // comma-separated genes (pasted text: "inlA, hly, qacH").
    for raw_line in ids_text.lines() {
        let line = raw_line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        // skip a header row on the first line
        if !any_entry && is_header_line(line) {
            continue;
        }
        for cell in line.split([',', ';']) {
            let cell = cell.trim().trim_matches('"').trim();
            if cell.is_empty() || cell.starts_with('#') {
                continue;
            }
            any_entry = true;
            let mut resolved: Option<(&Gene, String)> = None;
            // parenthesized locus tags first: "pva (lmo0446)"
            for m in parentheses(cell) {
                if let Some(g) = by_locus.get(m.as_str()) {
                    resolved = Some((g, g.locus_tag.clone()));
                    break;
                }
            }
            if resolved.is_none() {
                // then each whitespace-separated token of the cell
                for tok in cell.split_whitespace().filter(|t| !t.starts_with('#')) {
                    let bare = tok
                        .trim_start_matches(['(', '['])
                        .trim_end_matches([')', ']']);
                    if let Some(g) = by_locus.get(bare) {
                        resolved = Some((g, g.locus_tag.clone()));
                        break;
                    }
                    if let Some(g) = by_old_locus.get(bare) {
                        // re-annotated genome: keep the user's (old) spelling
                        resolved = Some((g, bare.to_string()));
                        break;
                    }
                    if let Some(gs) = by_symbol.get(&bare.to_lowercase()) {
                        let g = if gs.len() == 1 {
                            gs[0]
                        } else {
                            let cands = gs
                                .iter()
                                .map(|g| {
                                    Ok(Candidate {
                                        gene: g,
                                        seq: gene_sequence(ref_fasta, g)?,
                                    })
                                })
                                .collect::<Result<Vec<_>>>()?;
                            let (i, why) = choose(bare, &cands)
                                .filter(|(i, _)| *i < gs.len())
                                .unwrap_or((0, "the first in the genome".into()));
                            ambiguous.push(ambiguity_note(bare, gs, i, &why));
                            gs[i]
                        };
                        // keep the user's spelling as the FASTA header
                        resolved = Some((g, bare.to_string()));
                        break;
                    }
                }
            }
            match resolved {
                Some((gene, header)) => {
                    if used.contains(&header) {
                        continue; // same gene listed twice
                    }
                    used.push(header.clone());
                    found.push(header.clone());
                    let seq = gene_sequence(ref_fasta, gene)?;
                    fasta.push('>');
                    fasta.push_str(&header);
                    // where the sequence came from, for the run to show
                    fasta.push_str(" reference ");
                    fasta.push_str(&gene.locus_tag);
                    fasta.push('\n');
                    for chunk in seq.chunks(60) {
                        fasta.push_str(std::str::from_utf8(chunk).unwrap_or(""));
                        fasta.push('\n');
                    }
                }
                None => missing.push(cell.to_string()),
            }
        }
    }

    Ok(PanelFromIds {
        fasta,
        missing,
        found,
        ambiguous,
    })
}

/// "prfA names 2 genes in your reference: lmo0200 (listeriolysin positive
/// regulatory protein), taken (VFDB files it as prfA), and lmo2543
/// (peptide chain release factor 1). Write lmo2543 to search that one."
fn ambiguity_note(name: &str, genes: &[&Gene], chosen: usize, why: &str) -> String {
    let words = |g: &Gene| {
        if g.product.is_empty() {
            g.locus_tag.clone()
        } else {
            format!("{} ({})", g.locus_tag, g.product)
        }
    };
    let others: Vec<&&Gene> = genes
        .iter()
        .enumerate()
        .filter(|(i, _)| *i != chosen)
        .map(|(_, g)| g)
        .collect();
    format!(
        "{name} names {} genes in your reference: {}, taken ({why}), and {}. Write {} to search {}.",
        genes.len(),
        words(genes[chosen]),
        others.iter().map(|g| words(g)).collect::<Vec<_>>().join(" and "),
        others
            .iter()
            .map(|g| g.locus_tag.as_str())
            .collect::<Vec<_>>()
            .join(" or "),
        if others.len() == 1 { "that one" } else { "one of those" },
    )
}

fn is_header_line(line: &str) -> bool {
    let l = line.to_lowercase();
    [
        "gene",
        "gene_id",
        "genes",
        "id",
        "ids",
        "locus_tag",
        "locus",
        "name",
    ]
    .iter()
    .any(|h| l == *h || l.starts_with(&format!("{h},")))
}

/// Contents of all (...) groups in the line.
fn parentheses(line: &str) -> Vec<String> {
    let mut out = Vec::new();
    let bytes = line.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] == b'(' || bytes[i] == b'[' {
            let close = if bytes[i] == b'(' { b')' } else { b']' };
            let start = i + 1;
            let mut j = start;
            while j < bytes.len() && bytes[j] != close {
                j += 1;
            }
            if j > start {
                out.push(line[start..j].trim().to_string());
            }
            i = j + 1;
        } else {
            i += 1;
        }
    }
    out
}

fn gene_sequence(ref_fasta: &Path, gene: &Gene) -> Result<Vec<u8>> {
    // cached parse per call is fine: the panel is built once
    let recs = parse_fasta(ref_fasta)?;
    let rec = recs
        .iter()
        .find(|r| r.id == gene.seqid)
        .ok_or_else(|| {
            crate::friendly(format!(
                "The annotation mentions sequence \"{}\" which is not present in the reference genome file.",
                gene.seqid
            ))
        })?;
    let s = gene.start as usize - 1;
    let e = gene.end as usize;
    let seq = &rec.seq[s.min(rec.seq.len())..e.min(rec.seq.len())];
    if gene.strand < 0 {
        Ok(revcomp(seq))
    } else {
        Ok(seq.to_vec())
    }
}

/// A panel assembled from several sources, numbered again: each gene's
/// sequences become `gene`, `gene__v2`, ... in the order they came, and a
/// sequence a gene already holds is dropped. The reference copy, added
/// first, stays the first variant.
pub fn renumber_variants(records: &[crate::fasta::FastaRecord]) -> String {
    let mut genes: Vec<(&str, Vec<&crate::fasta::FastaRecord>)> = Vec::new();
    for r in records {
        let g = variant_gene(&r.id);
        match genes.iter_mut().find(|(name, _)| *name == g) {
            Some((_, held)) => {
                if held.iter().all(|h| !h.seq.eq_ignore_ascii_case(&r.seq)) {
                    held.push(r);
                }
            }
            None => genes.push((g, vec![r])),
        }
    }
    let mut out = String::new();
    for (g, held) in genes {
        for (i, r) in held.iter().enumerate() {
            out.push('>');
            out.push_str(&variant_id(g, i + 1));
            if !r.desc.is_empty() {
                out.push(' ');
                out.push_str(&r.desc);
            }
            out.push('\n');
            for chunk in r.seq.chunks(60) {
                out.push_str(&String::from_utf8_lossy(chunk));
                out.push('\n');
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// EGD-e's two "prfA": PrfA (lmo0200) first in the genome, RF1
    /// (lmo2543) later.
    fn two_prfa(tag: &str) -> (std::path::PathBuf, std::path::PathBuf) {
        let d = std::env::temp_dir().join(format!("sc-panel-{tag}-{}", std::process::id()));
        std::fs::create_dir_all(&d).unwrap();
        std::fs::write(d.join("r.fa"), ">c\nAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT\n").unwrap();
        std::fs::write(
            d.join("r.gff"),
            "##gff-version 3\n\
             c\tR\tgene\t1\t8\t.\t+\t.\tID=a;gene=prfA;locus_tag=lmo0200;gene_biotype=protein_coding\n\
             c\tR\tCDS\t1\t8\t.\t+\t0\tParent=a;locus_tag=lmo0200;product=listeriolysin positive regulatory protein\n\
             c\tR\tgene\t17\t24\t.\t+\t.\tID=b;gene=prfA;locus_tag=lmo2543;gene_biotype=protein_coding\n\
             c\tR\tCDS\t17\t24\t.\t+\t0\tParent=b;locus_tag=lmo2543;product=peptide chain release factor 1\n",
        )
        .unwrap();
        (d.join("r.fa"), d.join("r.gff"))
    }

    #[test]
    fn a_name_of_two_genes_takes_the_first_and_says_so() {
        let (fa, gff) = two_prfa("first");
        let p = panel_from_ids(&fa, &gff, "prfA").unwrap();
        // the first in genome order, not the last (which was RF1)
        assert!(
            p.fasta.starts_with(">prfA reference lmo0200\nAAAACCCC"),
            "{}",
            p.fasta
        );
        assert_eq!(p.ambiguous.len(), 1);
        let note = &p.ambiguous[0];
        assert!(note.contains("prfA names 2 genes"), "{note}");
        assert!(note.contains("lmo0200 (listeriolysin positive regulatory protein), taken (the first in the genome)"), "{note}");
        assert!(note.contains("Write lmo2543 to search that one."), "{note}");
        // asking by locus tag is never ambiguous
        let p = panel_from_ids(&fa, &gff, "lmo2543").unwrap();
        assert!(p.ambiguous.is_empty());
    }

    #[test]
    fn a_chooser_picks_among_genes_of_one_name() {
        let (fa, gff) = two_prfa("choose");
        let pick_rf1 = |name: &str, c: &[Candidate]| {
            assert_eq!(name, "prfA");
            assert_eq!(c[1].seq, b"AAAACCCC");
            c.iter()
                .position(|c| c.gene.product.contains("release factor"))
                .map(|i| (i, "a test said so".to_string()))
        };
        let p = panel_from_ids_with(&fa, &gff, "prfA", &pick_rf1).unwrap();
        assert!(
            p.fasta.starts_with(">prfA reference lmo2543"),
            "{}",
            p.fasta
        );
        assert!(p.ambiguous[0]
            .contains("lmo2543 (peptide chain release factor 1), taken (a test said so)"));
    }

    #[test]
    fn renumbers_variants_from_several_sources() {
        let recs = crate::fasta::parse_fasta_str(
            ">cadA reference lmo0001\nAAAA\n>inlA reference lmo0433\nCCCC\n\
             >cadA__v2 AMRFinderPlus cadA_Lm: CadA [X.1:1-4]\nGGGG\n\
             >cadA__v2 Library g11 cadA: ATPase [Listeria library t]\nTTTT\n\
             >cadA__v3 Library g10 cadA: same as the reference\naaaa\n",
        )
        .unwrap();
        let out = renumber_variants(&recs);
        let ids: Vec<&str> = out
            .lines()
            .filter_map(|l| l.strip_prefix('>'))
            .map(|l| l.split_whitespace().next().unwrap())
            .collect();
        assert_eq!(ids, ["cadA", "cadA__v2", "cadA__v3", "inlA"]);
        assert!(out.contains(">cadA__v3 Library g11 cadA: ATPase [Listeria library t]\nTTTT\n"));
    }

    #[test]
    fn variant_ids_name_their_gene() {
        assert_eq!(variant_gene("cadA"), "cadA");
        assert_eq!(variant_gene("cadA__v2"), "cadA");
        assert_eq!(variant_gene("cadA__vx"), "cadA__vx");
        assert_eq!(variant_gene("__v2"), "__v2");
        assert_eq!(variant_id("cadA", 1), "cadA");
        assert_eq!(variant_id("cadA", 3), "cadA__v3");
    }
}
