//! Unaligned gap regions on the reference and the genes inside them.

use crate::delta::DeltaFile;
use crate::gff::Gene;
use straincompass_types::GapRow;
use std::collections::HashMap;

/// Compute unaligned regions (>= min_gap bp) per reference sequence, with
/// the genes overlapping each region.
pub fn unaligned_gaps(
    genes: &[Gene],
    delta: &DeltaFile,
    ref_lengths: &HashMap<String, u64>,
    min_gap: u64,
) -> Vec<GapRow> {
    let aligned = delta.aligned_ref_intervals();
    let mut genes_by_seqid: HashMap<&str, Vec<&Gene>> = HashMap::new();
    for g in genes {
        genes_by_seqid.entry(g.seqid.as_str()).or_default().push(g);
    }

    let mut rows = Vec::new();
    let mut seqids: Vec<&String> = ref_lengths.keys().collect();
    seqids.sort();
    for seqid in seqids {
        let seqlen = ref_lengths[seqid];
        let mut ivs = aligned.get(seqid).cloned().unwrap_or_default();
        crate::delta::merge_intervals(&mut ivs);
        // complement
        let mut cursor = 1u64;
        let mut gaps: Vec<(u64, u64)> = Vec::new();
        for (s, e) in ivs {
            if s > cursor {
                gaps.push((cursor, s - 1));
            }
            cursor = e + 1;
        }
        if cursor <= seqlen {
            gaps.push((cursor, seqlen));
        }
        for (gs, ge) in gaps {
            let length = ge - gs + 1;
            if length < min_gap {
                continue;
            }
            let mut inside: Vec<&Gene> = Vec::new();
            if let Some(gs2) = genes_by_seqid.get(seqid.as_str()) {
                for g in gs2 {
                    if g.start <= ge && g.end >= gs {
                        inside.push(g);
                    }
                }
                inside.sort_by_key(|g| g.start);
            }
            rows.push(GapRow {
                seqid: seqid.clone(),
                start: gs,
                end: ge,
                length,
                n_genes: inside.len(),
                genes: inside.into_iter().map(|g| g.locus_tag.clone()).collect(),
            });
        }
    }
    rows.sort_by_key(|a| (a.seqid.clone(), a.start));
    rows
}
