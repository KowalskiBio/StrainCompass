//! Per-gene alignment coverage, presence calls and exact variant counts.

use crate::delta::{merge_intervals, Alignment, DeltaFile};
use crate::fasta::FastaRecord;
use crate::gff::Gene;
use crate::msa;
use crate::{friendly, Result};
use std::collections::HashMap;
use straincompass_types::{Call, GeneCoverageRow};

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
                (_, b'-') => (0, 1),       // insertion in reference
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

    let loci = query_loci(genes, delta);
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
        // R parity: the call is made on cov_pct rounded to one decimal
        // (a 94.96% gene is present in R's table), and with the default
        // partial_cov of 0 any aligned base keeps a gene out of ABSENT
        let cov_call = (cov_pct * 10.0).round() / 10.0;
        let call = if cov_call >= present_cov {
            Call::Present
        } else if cov_bp > 0 && cov_pct > partial_cov {
            Call::Partial
        } else {
            Call::Absent
        };
        rows.push(GeneCoverageRow {
            locus_tag: g.locus_tag.clone(),
            symbol: g.symbol.clone(),
            protein_id: g.protein_id.clone(),
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
            qry_loci: loci[idx].clone(),
        });
    }
    Ok(rows)
}

/// A run of paired columns in a delta block: reference positions
/// `r..r+len` align base for base to query positions from `q`, stepping
/// down on a reverse block.
struct PairedRun {
    r: u64,
    q: u64,
    len: u64,
}

/// The paired runs of one block, in reference order. Positions only, so
/// no sequence is needed.
fn paired_runs(a: &Alignment) -> Vec<PairedRun> {
    let adv = |q: u64, n: u64| {
        if a.qry_rev {
            q.saturating_sub(n)
        } else {
            q + n
        }
    };
    let mut runs = Vec::with_capacity(a.deltas.len() + 1);
    let mut r = a.ref_start;
    let mut q = if a.qry_rev { a.qry_hi } else { a.qry_lo };
    for &d in &a.deltas {
        let k = d.unsigned_abs().saturating_sub(1);
        if k > 0 {
            runs.push(PairedRun { r, q, len: k });
            r += k;
            q = adv(q, k);
        }
        if d < 0 {
            q = adv(q, 1); // query base facing a reference gap
        } else {
            r += 1; // reference base facing a query gap
        }
    }
    if r <= a.ref_end {
        runs.push(PairedRun {
            r,
            q,
            len: a.ref_end - r + 1,
        });
    }
    runs
}

/// Where each gene lies in the query, aligned with `genes`: for every
/// query contig and strand its aligned bases map to, the query span as
/// "contig:start-end(+|-)", joined with "; " with the most aligned bases
/// first. A gene split over contigs (or duplicated) lists each place;
/// blocks landing next to each other on the same contig and strand are
/// one place, and slivers under a tenth of the best place are dropped.
/// Empty for a gene with nothing aligned.
pub fn query_loci(genes: &[Gene], delta: &DeltaFile) -> Vec<String> {
    let mut genes_by_seqid: HashMap<&str, Vec<usize>> = HashMap::new();
    for (idx, g) in genes.iter().enumerate() {
        genes_by_seqid
            .entry(g.seqid.as_str())
            .or_default()
            .push(idx);
    }
    struct Place<'a> {
        contig: &'a str,
        rev: bool,
        lo: u64,
        hi: u64,
        bp: u64,
    }
    let mut places: HashMap<usize, Vec<Place>> = HashMap::new();
    for a in &delta.alignments {
        let Some(gidx) = genes_by_seqid.get(a.ref_seqid.as_str()) else {
            continue;
        };
        let mut runs: Option<Vec<PairedRun>> = None;
        for &gi in gidx {
            let g = &genes[gi];
            if g.start > a.ref_end || g.end < a.ref_start {
                continue;
            }
            let runs = runs.get_or_insert_with(|| paired_runs(a));
            let first = runs.partition_point(|run| run.r + run.len <= g.start);
            let (mut lo, mut hi, mut bp) = (u64::MAX, 0u64, 0u64);
            for run in &runs[first..] {
                if run.r > g.end {
                    break;
                }
                let ra = run.r.max(g.start);
                let rb = (run.r + run.len - 1).min(g.end);
                let at = |p: u64| {
                    if a.qry_rev {
                        run.q.saturating_sub(p - run.r)
                    } else {
                        run.q + (p - run.r)
                    }
                };
                let (qa, qb) = (at(ra), at(rb));
                lo = lo.min(qa.min(qb));
                hi = hi.max(qa.max(qb));
                bp += rb - ra + 1;
            }
            if bp == 0 {
                continue;
            }
            let list = places.entry(gi).or_default();
            let slack = g.end - g.start + 1;
            match list.iter_mut().find(|p| {
                p.contig == a.qry_seqid
                    && p.rev == a.qry_rev
                    && lo <= p.hi + slack
                    && hi + slack >= p.lo
            }) {
                Some(p) => {
                    p.lo = p.lo.min(lo);
                    p.hi = p.hi.max(hi);
                    p.bp += bp;
                }
                None => list.push(Place {
                    contig: &a.qry_seqid,
                    rev: a.qry_rev,
                    lo,
                    hi,
                    bp,
                }),
            }
        }
    }
    let mut out = vec![String::new(); genes.len()];
    for (gi, mut list) in places {
        list.sort_by(|x, y| {
            y.bp.cmp(&x.bp)
                .then(x.contig.cmp(y.contig))
                .then(x.lo.cmp(&y.lo))
        });
        // a block edge spilling a few bases into the gene is not a place
        // of it: keep the places carrying at least a tenth of the best
        let floor = list[0].bp / 10;
        list.retain(|p| p.bp >= floor);
        out[gi] = list
            .iter()
            .map(|p| {
                let strand = if p.rev { '-' } else { '+' };
                format!("{}:{}-{}({strand})", p.contig, p.lo, p.hi)
            })
            .collect::<Vec<_>>()
            .join("; ");
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::delta::parse_delta_str;

    fn gene(tag: &str, start: u64, end: u64) -> Gene {
        Gene {
            seqid: "ref".into(),
            start,
            end,
            strand: 1,
            locus_tag: tag.into(),
            old_locus_tag: String::new(),
            symbol: String::new(),
            biotype: String::new(),
            protein_id: String::new(),
            product: String::new(),
        }
    }

    #[test]
    fn maps_genes_through_forward_and_reverse_blocks() {
        // ref 1-100 on ctgA 1001-1100 forward; ref 201-300 on ctgB
        // 500-401 reverse, with an extra query base after ref 210 and
        // ref 250 facing a gap in the query
        let delta = parse_delta_str(
            "/r.fa /q.fa\nNUCMER\n>ref ctgA 1000 5000\n1 100 1001 1100 0 0 0\n0\n\
             >ref ctgB 1000 1000\n201 300 500 401 2 2 0\n-11\n40\n0\n",
        )
        .unwrap();
        let genes = vec![
            gene("fwd", 11, 30),
            gene("rev", 221, 230),
            gene("none", 150, 180),
        ];
        let loci = query_loci(&genes, &delta);
        assert_eq!(loci[0], "ctgA:1011-1030(+)");
        // ref 221 sits 11 columns past the insertion: query 500 - 10 - 1
        // (the inserted base) - 10 = 479, and ref 230 is 470
        assert_eq!(loci[1], "ctgB:470-479(-)");
        assert_eq!(loci[2], "");
    }
}
