//! Per-gene alignment coverage, presence calls and exact variant counts.

use crate::delta::{merge_intervals, DeltaFile};
use crate::fasta::FastaRecord;
use crate::gff::Gene;
use crate::msa;
use crate::{friendly, Result};
use bactiment_types::{Call, GeneCoverageRow};
use std::collections::HashMap;

/// Compute the genes coverage table for one query.
///
/// `ref_records` are the sanitized reference sequences; `qry_records` the
/// sanitized query sequences.
pub fn gene_coverage(
    genes: &[Gene],
    delta: &DeltaFile,
    ref_records: &[FastaRecord],
    qry_records: &[FastaRecord],
    present_cov: f64,
    partial_cov: f64,
) -> Result<Vec<GeneCoverageRow>> {
    let ref_by_id: HashMap<&str, &FastaRecord> =
        ref_records.iter().map(|r| (r.id.as_str(), r)).collect();
    let qry_by_id: HashMap<&str, &FastaRecord> =
        qry_records.iter().map(|r| (r.id.as_str(), r)).collect();

    // Per gene accumulation.
    struct Acc {
        cov_bp: u64,
        best_identity: f64,
        mismatches: u64,
        indels: u64,
    }
    let mut acc: HashMap<usize, Acc> = HashMap::new();

    // Coverage from merged aligned intervals, per reference seqid.
    let aligned = delta.aligned_ref_intervals();

    // Index genes per seqid for overlap search.
    let mut genes_by_seqid: HashMap<&str, Vec<usize>> = HashMap::new();
    for (idx, g) in genes.iter().enumerate() {
        genes_by_seqid
            .entry(g.seqid.as_str())
            .or_default()
            .push(idx);
    }

    for (seqid, ivs) in &aligned {
        let Some(gidx) = genes_by_seqid.get(seqid.as_str()) else {
            continue;
        };
        let mut ivs = ivs.clone();
        merge_intervals(&mut ivs);
        // Walk genes and intervals together (both sorted by start).
        let mut gene_list: Vec<(u64, u64, usize)> = gidx
            .iter()
            .map(|&i| (genes[i].start, genes[i].end, i))
            .collect();
        gene_list.sort();
        for &(gstart, gend, gi) in &gene_list {
            let mut cov = 0u64;
            for &(s, e) in &ivs {
                if e < gstart {
                    continue;
                }
                if s > gend {
                    break;
                }
                cov += e.min(gend).saturating_sub(s.max(gstart)) + 1;
            }
            if cov > 0 {
                acc.insert(
                    gi,
                    Acc {
                        cov_bp: cov,
                        best_identity: 0.0,
                        mismatches: 0,
                        indels: 0,
                    },
                );
            }
        }
    }

    // Exact variant stats + best identity from block reconstruction.
    for a in &delta.alignments {
        let Some(gidx) = genes_by_seqid.get(a.ref_seqid.as_str()) else {
            continue;
        };
        // Only blocks that overlap at least one gene matter.
        let overlapping: Vec<usize> = gidx
            .iter()
            .copied()
            .filter(|&i| {
                let g = &genes[i];
                g.start <= a.ref_end && g.end >= a.ref_start
            })
            .collect();
        if overlapping.is_empty() {
            continue;
        }
        let Some(ref_rec) = ref_by_id.get(a.ref_seqid.as_str()) else {
            return Err(friendly(format!(
                "The reference sequence {} was not found in the reference FASTA file.",
                a.ref_seqid
            )));
        };
        let Some(qry_rec) = qry_by_id.get(a.qry_seqid.as_str()) else {
            continue;
        };
        let (ref_bases, qry_bases) = msa::block_bases(a, ref_rec, qry_rec);
        let pw = msa::reconstruct(a, &ref_bases, &qry_bases);
        let idy = a.identity();
        // For each column, attribute to genes containing the ref position.
        // Collect events per gene: (mismatches, indels).
        let mut col_events: HashMap<usize, (u64, u64)> = HashMap::new();
        for c in 0..pw.len() {
            let rp = pw.ref_pos[c];
            if rp < a.ref_start {
                continue;
            }
            let (m, ind) = match (pw.ref_row[c], pw.qry_row[c]) {
                (b'-', b'-') => continue,
                (b'-', _) => (0u64, 1u64), // insertion in query
                (_, b'-') => (0, 1),      // insertion in reference
                (r, q) => {
                    if r != q {
                        (1, 0)
                    } else {
                        continue;
                    }
                }
            };
            for &gi in &overlapping {
                let g = &genes[gi];
                if rp >= g.start && rp <= g.end {
                    let e = col_events.entry(gi).or_insert((0, 0));
                    e.0 += m;
                    e.1 += ind;
                    break; // a position belongs to one gene (non-overlapping assumed)
                }
            }
        }
        for (gi, (m, ind)) in col_events {
            let e = acc.entry(gi).or_insert(Acc {
                cov_bp: 0,
                best_identity: 0.0,
                mismatches: 0,
                indels: 0,
            });
            e.mismatches += m;
            e.indels += ind;
            e.best_identity = e.best_identity.max(idy);
        }
        // Best identity also for genes that are fully covered without variants.
        for &gi in &overlapping {
            if let Some(e) = acc.get_mut(&gi) {
                e.best_identity = e.best_identity.max(idy);
            }
        }
    }

    let mut rows = Vec::with_capacity(genes.len());
    for (idx, g) in genes.iter().enumerate() {
        let length = g.end - g.start + 1;
        let (cov_bp, best_identity, mismatches, indels) = match acc.get(&idx) {
            Some(a) => (a.cov_bp, a.best_identity, a.mismatches, a.indels),
            None => (0, 0.0, 0, 0),
        };
        let cov_pct = if length > 0 {
            100.0 * cov_bp as f64 / length as f64
        } else {
            0.0
        };
        let call = if cov_pct >= present_cov {
            Call::Present
        } else if cov_pct > partial_cov {
            Call::Partial
        } else {
            Call::Absent
        };
        rows.push(GeneCoverageRow {
            locus_tag: g.locus_tag.clone(),
            symbol: g.symbol.clone(),
            biotype: g.biotype.clone(),
            seqid: g.seqid.clone(),
            start: g.start,
            end: g.end,
            length,
            cov_bp,
            cov_pct,
            call,
            best_identity,
            mismatches,
            indels,
        });
    }
    Ok(rows)
}
