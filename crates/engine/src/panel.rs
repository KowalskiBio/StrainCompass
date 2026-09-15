//! Build a gene panel FASTA from a list of gene identifiers.
//!
//! Mirrors how `genes_of_interest.fasta` was produced for the R pipeline:
//! the user supplies identifiers (locus tags like lmo0444, gene symbols
//! like inlA, or "symbol (locus_tag)" pairs), and each gene's sequence is
//! extracted from the reference genome in the gene's own orientation.

use crate::fasta::{parse_fasta, revcomp};
use crate::gff::{parse_gff, Gene};
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
}

/// Extract gene sequences from the reference.
///
/// `ids_text` is the raw content of the user's CSV/TSV/list file; every
/// line holds one gene (extra columns and comments are tolerated, and the
/// parenthesized locus tag in entries like "pva (lmo0446)" wins over the
/// symbol).
pub fn panel_from_ids(ref_fasta: &Path, ref_gff: &Path, ids_text: &str) -> Result<PanelFromIds> {
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
    let mut by_symbol: HashMap<String, &Gene> = HashMap::new();
    for g in &genes {
        by_locus.insert(g.locus_tag.as_str(), g);
        if !g.old_locus_tag.is_empty() {
            by_old_locus.insert(g.old_locus_tag.as_str(), g);
        }
        if !g.symbol.is_empty() {
            by_symbol.insert(g.symbol.to_lowercase(), g);
        }
    }

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
                    if let Some(g) = by_symbol.get(&bare.to_lowercase()) {
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
    })
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
