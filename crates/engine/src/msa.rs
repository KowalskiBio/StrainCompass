//! Pairwise alignment reconstruction from the delta file, used by the
//! gene alignment (MSA) viewer and for exact per-gene variant counts.

use crate::delta::Alignment;
use crate::fasta::{revcomp, FastaRecord};

/// A reconstructed pairwise alignment of one delta block.
pub struct Pairwise {
    /// Reference row (bases or '-'), forward reference orientation.
    pub ref_row: Vec<u8>,
    /// Query row (bases or '-'), forward query orientation.
    pub qry_row: Vec<u8>,
    /// Reference position (1-based) of the last emitted reference base at
    /// each column; for gap columns in the reference row this repeats the
    /// previous position (ref_start - 1 before the first base).
    pub ref_pos: Vec<u64>,
    /// Query position (1-based, forward orientation) of the last emitted
    /// query base at each column (qry_lo - 1 before the first base).
    pub qry_pos: Vec<u64>,
}

impl Pairwise {
    pub fn len(&self) -> usize {
        self.ref_row.len()
    }
    pub fn is_empty(&self) -> bool {
        self.ref_row.is_empty()
    }
}

/// Reconstruct the full pairwise alignment of one delta block.
///
/// `ref_bases` must be reference[res_start-1 .. res_end] (forward);
/// `qry_bases` must be the query bases qry_lo..=qry_hi, reverse
/// complemented when the block is on the reverse strand.
pub fn reconstruct(a: &Alignment, ref_bases: &[u8], qry_bases: &[u8]) -> Pairwise {
    let cap = a.aligned_len() as usize;
    let mut pw = Pairwise {
        ref_row: Vec::with_capacity(cap),
        qry_row: Vec::with_capacity(cap),
        ref_pos: Vec::with_capacity(cap),
        qry_pos: Vec::with_capacity(cap),
    };

    struct Walk<'a> {
        a: &'a Alignment,
        ref_bases: &'a [u8],
        qry_bases: &'a [u8],
        i: usize,
        j: usize,
        r: u64,
        q: u64,
    }

    impl Walk<'_> {
        fn next_q(&mut self) -> u64 {
            if self.a.qry_rev {
                self.q -= 1;
            } else {
                self.q += 1;
            }
            self.q
        }
        fn emit_paired(&mut self, k: usize, pw: &mut Pairwise) {
            for _ in 0..k {
                if self.i >= self.ref_bases.len() || self.j >= self.qry_bases.len() {
                    return;
                }
                self.r += 1;
                let q = self.next_q();
                pw.ref_row.push(self.ref_bases[self.i]);
                pw.qry_row.push(self.qry_bases[self.j]);
                pw.ref_pos.push(self.r);
                pw.qry_pos.push(q);
                self.i += 1;
                self.j += 1;
            }
        }
    }

    let mut w = Walk {
        a,
        ref_bases,
        qry_bases,
        i: 0,
        j: 0,
        r: a.ref_start - 1,
        q: if a.qry_rev {
            a.qry_hi + 1
        } else {
            a.qry_lo - 1
        },
    };

    for &d in &a.deltas {
        let k = d.unsigned_abs() as usize;
        w.emit_paired(k.saturating_sub(1), &mut pw);
        if d < 0 {
            // insertion in the query: gap in the reference row
            if w.j < qry_bases.len() {
                let q = w.next_q();
                pw.ref_row.push(b'-');
                pw.qry_row.push(qry_bases[w.j]);
                pw.ref_pos.push(w.r);
                pw.qry_pos.push(q);
                w.j += 1;
            }
        } else {
            // insertion in the reference: gap in the query row
            if w.i < ref_bases.len() {
                w.r += 1;
                pw.ref_row.push(ref_bases[w.i]);
                pw.qry_row.push(b'-');
                pw.ref_pos.push(w.r);
                pw.qry_pos.push(w.q);
                w.i += 1;
            }
        }
    }
    // remaining paired bases
    w.emit_paired(usize::MAX, &mut pw);

    pw
}

/// Extract query bases covering a reference interval, from a block.
/// Returns (qry_start, qry_end, aligned_columns) where coordinates are in
/// the query's forward orientation, or None if the block does not cover
/// the interval.
pub struct GeneSlice {
    pub ref_row: Vec<u8>,
    pub qry_row: Vec<u8>,
    /// 1-based inclusive reference span covered.
    pub ref_start: u64,
    pub ref_end: u64,
    /// Query span covered (forward orientation).
    pub qry_start: u64,
    pub qry_end: u64,
}

pub fn slice_block_to_ref_range(
    pw: &Pairwise,
    a: &Alignment,
    ref_from: u64,
    ref_to: u64,
) -> GeneSlice {
    let mut ref_row = Vec::new();
    let mut qry_row = Vec::new();
    let mut first_q: Option<u64> = None;
    let mut last_q: u64 = 0;
    let mut ref_pos_of_first: u64 = 0;
    let mut ref_pos_of_last: u64 = 0;
    for c in 0..pw.len() {
        let rp = pw.ref_pos[c];
        // A column belongs to the range if its reference position (last
        // emitted base) is inside; ref gap columns inherit the previous
        // position.
        if rp >= ref_from && rp <= ref_to {
            if ref_row.is_empty() {
                ref_pos_of_first = rp;
            }
            ref_pos_of_last = rp;
            ref_row.push(pw.ref_row[c]);
            qry_row.push(pw.qry_row[c]);
            if pw.qry_row[c] != b'-' {
                let qp = pw.qry_pos[c];
                first_q = Some(first_q.map_or(qp, |f: u64| f.min(qp)));
                last_q = last_q.max(qp);
            }
        }
    }
    let qry_start = first_q.unwrap_or(0);
    let qry_end = last_q;
    let _ = a;
    // The span actually covered, not the range that was asked for: a gene
    // whose alignment starts partway in would otherwise claim the whole gene
    // in its block header while carrying only the aligned columns, and the
    // unaligned stretches beside it would overlap it instead of tiling it.
    // Both stay 0 on an empty slice, which the caller skips.
    let (ref_start, ref_end) = if ref_row.is_empty() {
        (0, 0)
    } else {
        (ref_pos_of_first, ref_pos_of_last)
    };
    GeneSlice {
        ref_row,
        qry_row,
        ref_start,
        ref_end,
        qry_start,
        qry_end,
    }
}

/// Concatenate per-column variant stats over a reference range.
pub struct VariantStats {
    pub mismatches: u64,
    /// Bases in the query not present in the reference.
    pub qry_ins: u64,
    /// Bases in the reference not present in the query.
    pub ref_ins: u64,
}

pub fn variant_stats(pw: &Pairwise, ref_from: u64, ref_to: u64) -> VariantStats {
    let mut st = VariantStats {
        mismatches: 0,
        qry_ins: 0,
        ref_ins: 0,
    };
    for c in 0..pw.len() {
        let rp = pw.ref_pos[c];
        if rp < ref_from || rp > ref_to {
            continue;
        }
        match (pw.ref_row[c], pw.qry_row[c]) {
            (b'-', b'-') => {}
            (b'-', _) => st.qry_ins += 1,
            (_, b'-') => st.ref_ins += 1,
            (r, q) => {
                if r != q {
                    st.mismatches += 1;
                }
            }
        }
    }
    st
}

/// Translate a nucleotide sequence (already in gene orientation, i.e.
/// 5'->3' of the gene) and list premature stop codons.
pub fn premature_stops(seq: &[u8]) -> Vec<(u64, u64)> {
    let mut stops = Vec::new();
    let n = seq.len() / 3;
    for ci in 0..n {
        let codon = &seq[ci * 3..ci * 3 + 3];
        if is_stop(codon) && ci + 1 < n {
            // a stop before the last complete codon is premature
            stops.push((ci as u64, ci as u64 + 1));
        }
    }
    stops
}

fn is_stop(codon: &[u8]) -> bool {
    matches!(
        (
            codon[0].to_ascii_uppercase(),
            codon[1].to_ascii_uppercase(),
            codon[2].to_ascii_uppercase()
        ),
        (b'T', b'A', b'A') | (b'T', b'A', b'G') | (b'T', b'G', b'A')
    )
}

/// Orient a slice for a gene on the minus strand: reverse complement both
/// rows (gaps stay gaps).
pub fn orient_for_strand(strand: i8, ref_row: Vec<u8>, qry_row: Vec<u8>) -> (Vec<u8>, Vec<u8>) {
    if strand >= 0 {
        return (ref_row, qry_row);
    }
    (revcomp_gaps(&ref_row), revcomp_gaps(&qry_row))
}

fn revcomp_gaps(row: &[u8]) -> Vec<u8> {
    row.iter()
        .rev()
        .map(|&b| if b == b'-' { b'-' } else { revcomp(&[b])[0] })
        .collect()
}

/// Convenience: fetch block bases for reconstruction.
pub fn block_bases<'a>(
    a: &Alignment,
    ref_rec: &'a FastaRecord,
    qry_rec: &'a FastaRecord,
) -> (Vec<u8>, Vec<u8>) {
    let ref_bases = crate::fasta::subseq(ref_rec, a.ref_start, a.ref_end, false);
    let qry_bases = crate::fasta::subseq(qry_rec, a.qry_lo, a.qry_hi, a.qry_rev);
    (ref_bases, qry_bases)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn alignment(ref_start: u64, ref_end: u64, qry_rev: bool) -> Alignment {
        Alignment {
            ref_seqid: "ctg".into(),
            qry_seqid: "q".into(),
            ref_start,
            ref_end,
            qry_lo: 1,
            qry_hi: ref_end - ref_start + 1,
            qry_rev,
            errors: 0,
            deltas: vec![],
        }
    }

    /// A gene overlapping the alignment only partly must get the span that
    /// is actually covered, not the gene's own range: the block used to
    /// claim the whole gene while carrying only the aligned columns, which
    /// made it overlap the unaligned stretches instead of tiling them.
    #[test]
    fn slice_reports_the_covered_span_not_the_requested_one() {
        let a = alignment(10, 25, false);
        let ref_bases = b"AAAAAAAAAAAAAAAA".to_vec();
        let qry_bases = b"CCCCCCCCCCCCCCCC".to_vec();
        let pw = reconstruct(&a, &ref_bases, &qry_bases);

        // gene 5..30: alignment covers 10..25 of it
        let s = slice_block_to_ref_range(&pw, &a, 5, 30);
        assert_eq!(s.ref_row.len(), 16);
        assert_eq!(
            s.ref_start, 10,
            "span must start at the alignment, not the gene"
        );
        assert_eq!(
            s.ref_end, 25,
            "span must end at the alignment, not the gene"
        );

        // gene 15..20: fully inside the alignment
        let s = slice_block_to_ref_range(&pw, &a, 15, 20);
        assert_eq!(s.ref_row.len(), 6);
        assert_eq!((s.ref_start, s.ref_end), (15, 20));

        // no overlap at all: empty slice, caller skips it
        let s = slice_block_to_ref_range(&pw, &a, 40, 50);
        assert!(s.ref_row.is_empty());
    }
}
