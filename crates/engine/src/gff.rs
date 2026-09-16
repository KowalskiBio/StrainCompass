//! Minimal GFF3 reader for reference annotations.

use crate::{friendly, Result};
use std::collections::HashMap;
use std::io::Read;
use std::path::Path;

#[derive(Debug, Clone)]
pub struct Gene {
    pub seqid: String,
    /// 1-based inclusive.
    pub start: u64,
    pub end: u64,
    pub strand: i8,
    pub locus_tag: String,
    /// Older tag kept in `old_locus_tag` by re-annotated genomes.
    pub old_locus_tag: String,
    pub symbol: String,
    pub biotype: String,
    /// RefSeq/GenBank protein accession (`protein_id=WP_...`) from the CDS
    /// child, joined in by locus tag. Empty when the gene has no coding child.
    pub protein_id: String,
    /// Function annotation from the `product` attribute (of the gene
    /// feature itself, or of the CDS feature sharing its locus tag).
    pub product: String,
}

/// Parse a GFF3 file into a list of genes.
///
/// Prefers `gene` / `pseudogene` features (NCBI style). If the file has no
/// gene-level features, falls back to merging CDS / tRNA / rRNA / ncRNA
/// features that share a locus tag.
pub fn parse_gff<P: AsRef<Path>>(path: P) -> Result<Vec<Gene>> {
    let mut text = String::new();
    std::fs::File::open(path.as_ref())?.read_to_string(&mut text)?;
    parse_gff_str(&text)
}

fn parse_attrs(field: &str) -> HashMap<String, String> {
    let mut map = HashMap::new();
    for kv in field.split(';') {
        let kv = kv.trim();
        if kv.is_empty() {
            continue;
        }
        if let Some(eq) = kv.find('=') {
            let (k, v) = kv.split_at(eq);
            map.insert(k.trim().to_string(), v[1..].to_string());
        }
    }
    map
}

pub fn parse_gff_str(text: &str) -> Result<Vec<Gene>> {
    let mut gene_rows: Vec<Gene> = Vec::new();
    // (locus tag, seqid, start, end, strand, feature type, attrs)
    type CdsRow = (
        String,
        String,
        u64,
        u64,
        i8,
        String,
        HashMap<String, String>,
    );
    let mut cds_rows: Vec<CdsRow> = Vec::new();
    let mut saw_gff_header = false;
    let mut n_data_rows = 0;

    for (lineno, line) in text.lines().enumerate() {
        let line = line.trim_end_matches('\r');
        if line.starts_with("##") {
            if line.starts_with("##gff-version") {
                saw_gff_header = true;
            }
            continue;
        }
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 9 {
            if n_data_rows == 0 && !saw_gff_header {
                return Err(friendly(
                    "This does not look like a GFF3 file. A GFF file has 9 columns separated by tabs and starts with \u{201c}##gff-version 3\u{201d}.".to_string()
                ));
            }
            continue;
        }
        n_data_rows += 1;
        let seqid = cols[0];
        let ftype = cols[2];
        let (start, end): (u64, u64) = match (cols[3].parse(), cols[4].parse()) {
            (Ok(s), Ok(e)) => (s, e),
            _ => {
                return Err(friendly(format!(
                    "Line {} of the GFF file has coordinates that are not numbers.",
                    lineno + 1
                )));
            }
        };
        let strand = match cols[6] {
            "+" => 1,
            "-" => -1,
            _ => 0,
        };
        let attrs = parse_attrs(cols[8]);

        match ftype {
            // feature type "gene" only, same as the R pipeline: GFFs mark
            // pseudogenes with their own feature type and R excluded them
            "gene" => {
                let locus_tag = attrs
                    .get("locus_tag")
                    .cloned()
                    .unwrap_or_else(|| format!("{}_{}_{}_gene", seqid, start, end));
                // R parity: symbol comes from the gene= attribute only
                let symbol = attrs.get("gene").cloned().unwrap_or_default();
                let biotype = attrs
                    .get("gene_biotype")
                    .cloned()
                    .unwrap_or_else(|| ftype.to_string());
                gene_rows.push(Gene {
                    seqid: seqid.to_string(),
                    start,
                    end,
                    strand,
                    locus_tag,
                    old_locus_tag: attrs.get("old_locus_tag").cloned().unwrap_or_default(),
                    symbol,
                    biotype,
                    protein_id: String::new(),
                    product: attrs.get("product").cloned().unwrap_or_default(),
                });
            }
            "CDS" | "tRNA" | "rRNA" | "ncRNA" | "tmRNA" => {
                cds_rows.push((
                    attrs.get("locus_tag").cloned().unwrap_or_default(),
                    seqid.to_string(),
                    start,
                    end,
                    strand,
                    ftype.to_string(),
                    attrs,
                ));
            }
            _ => {}
        }
    }

    if n_data_rows == 0 {
        return Err(friendly(
            "No annotation features were found in this GFF file.",
        ));
    }

    if !gene_rows.is_empty() {
        // Gene features usually carry no product; join it in from the CDS
        // features that share their locus tag. The protein accession is on
        // the CDS child too, same join.
        let mut products: HashMap<String, String> = HashMap::new();
        let mut proteins: HashMap<String, String> = HashMap::new();
        for (tag, _, _, _, _, _, attrs) in &cds_rows {
            if tag.is_empty() {
                continue;
            }
            if let Some(p) = attrs.get("product") {
                if !p.is_empty() {
                    products.entry(tag.clone()).or_insert_with(|| p.clone());
                }
            }
            if let Some(p) = attrs.get("protein_id") {
                if !p.is_empty() {
                    proteins.entry(tag.clone()).or_insert_with(|| p.clone());
                }
            }
        }
        for g in &mut gene_rows {
            if g.product.is_empty() {
                g.product = products.get(&g.locus_tag).cloned().unwrap_or_default();
            }
            if g.protein_id.is_empty() {
                g.protein_id = proteins.get(&g.locus_tag).cloned().unwrap_or_default();
            }
        }
        gene_rows.sort_by(|a, b| {
            (a.seqid.clone(), a.start, a.end).cmp(&(b.seqid.clone(), b.start, b.end))
        });
        return Ok(gene_rows);
    }

    // Fallback: merge same-locus_tag features into gene spans.
    let mut by_tag: HashMap<String, Gene> = HashMap::new();
    let mut order: Vec<String> = Vec::new();
    for (tag, seqid, start, end, strand, ftype, attrs) in cds_rows {
        let locus_tag = if !tag.is_empty() {
            tag
        } else {
            format!("{}_{}_{}_{}", seqid, start, end, "feat")
        };
        let entry = by_tag.entry(locus_tag.clone()).or_insert_with(|| {
            order.push(locus_tag.clone());
            Gene {
                seqid: seqid.clone(),
                start,
                end,
                strand,
                locus_tag: locus_tag.clone(),
                old_locus_tag: attrs.get("old_locus_tag").cloned().unwrap_or_default(),
                symbol: attrs
                    .get("gene")
                    .or_else(|| attrs.get("Name"))
                    .cloned()
                    .unwrap_or_default(),
                biotype: ftype.clone(),
                protein_id: attrs.get("protein_id").cloned().unwrap_or_default(),
                product: attrs.get("product").cloned().unwrap_or_default(),
            }
        });
        entry.start = entry.start.min(start);
        entry.end = entry.end.max(end);
        // Prefer a richer biotype when seen later (e.g. tRNA over CDS parts).
        if entry.biotype == "CDS" && ftype != "CDS" {
            entry.biotype = ftype.clone();
        }
        if entry.product.is_empty() {
            if let Some(p) = attrs.get("product") {
                if !p.is_empty() {
                    entry.product = p.clone();
                }
            }
        }
        if entry.protein_id.is_empty() {
            if let Some(p) = attrs.get("protein_id") {
                if !p.is_empty() {
                    entry.protein_id = p.clone();
                }
            }
        }
    }
    let mut genes: Vec<Gene> = order
        .into_iter()
        .filter_map(|t| by_tag.remove(&t))
        .collect();
    genes.sort_by_key(|a| (a.seqid.clone(), a.start, a.end));
    Ok(genes)
}
