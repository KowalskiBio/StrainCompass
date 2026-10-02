//! Base-level alignment of a panel gene to where it was found in a strain,
//! for the side-by-side "why does it differ" view: a Present strain next
//! to a Partial one shows exactly where the two copies part ways.
//!
//! Local alignment (Smith-Waterman, affine gaps) with blastn's default
//! scores, so the picture agrees with the search that made the call: the
//! aligned stretch is what BLAST would report, and the parts of the gene
//! outside it are what the strain does not share.

use crate::fasta::revcomp;

/// blastn defaults (-task blastn): +2/-3, gap open 5, extend 2.
const MATCH: i32 = 2;
const MISMATCH: i32 = -3;
const GAP_OPEN: i32 = 5;
const GAP_EXTEND: i32 = 2;

/// The alignment of a panel gene against a stretch of a strain contig.
#[derive(Debug, Clone, PartialEq)]
pub struct PanelAlignment {
    /// Aligned rows, '-' for gaps; panel gene on top.
    pub panel_row: Vec<u8>,
    pub strain_row: Vec<u8>,
    /// 1-based inclusive range of the panel gene that aligned.
    pub panel_start: u64,
    pub panel_end: u64,
    /// 1-based inclusive range of the window that aligned (window
    /// orientation, i.e. already reverse-complemented for a minus hit).
    pub window_start: u64,
    pub window_end: u64,
    pub identity: f64,
    pub mismatches: u64,
    pub gap_bases: u64,
}

impl PanelAlignment {
    /// Share of the panel gene inside the aligned stretch, percent.
    pub fn panel_coverage(&self, panel_len: usize) -> f64 {
        if panel_len == 0 || self.panel_end < self.panel_start {
            return 0.0;
        }
        (self.panel_end - self.panel_start + 1) as f64 * 100.0 / panel_len as f64
    }
}

#[derive(Clone, Copy, PartialEq)]
enum From {
    Stop,
    Diag,
    Up,
    Left,
}

/// Best local alignment of `panel` to `window`. None when nothing scores
/// above zero (no shared stretch at all).
pub fn align(panel: &[u8], window: &[u8]) -> Option<PanelAlignment> {
    let (n, m) = (panel.len(), window.len());
    if n == 0 || m == 0 {
        return None;
    }
    let w = m + 1;
    let neg = i32::MIN / 4;
    // H: best ending here; E: ending in a gap in the panel (Left), F: in
    // the window (Up). Traceback kept per matrix.
    let mut h_prev = vec![0i32; w];
    let mut h_cur = vec![0i32; w];
    let mut e_cur = vec![neg; w];
    let mut f_prev = vec![neg; w];
    let mut f_cur = vec![neg; w];
    let mut tb_h = vec![From::Stop; (n + 1) * w];
    // whether the gap in E/F at a cell continued an earlier gap
    let mut tb_e_ext = vec![false; (n + 1) * w];
    let mut tb_f_ext = vec![false; (n + 1) * w];
    let (mut best, mut bi, mut bj) = (0i32, 0usize, 0usize);
    for i in 1..=n {
        e_cur[0] = neg;
        h_cur[0] = 0;
        let a = panel[i - 1].to_ascii_uppercase();
        for j in 1..=m {
            let open_e = h_cur[j - 1] - GAP_OPEN - GAP_EXTEND;
            let ext_e = e_cur[j - 1] - GAP_EXTEND;
            e_cur[j] = open_e.max(ext_e);
            tb_e_ext[i * w + j] = ext_e > open_e;
            let open_f = h_prev[j] - GAP_OPEN - GAP_EXTEND;
            let ext_f = f_prev[j] - GAP_EXTEND;
            f_cur[j] = open_f.max(ext_f);
            tb_f_ext[i * w + j] = ext_f > open_f;
            let b = window[j - 1].to_ascii_uppercase();
            let s = if a == b && a != b'N' { MATCH } else { MISMATCH };
            let diag = h_prev[j - 1] + s;
            let (mut v, mut from) = (0, From::Stop);
            if diag > v {
                v = diag;
                from = From::Diag;
            }
            if f_cur[j] > v {
                v = f_cur[j];
                from = From::Up;
            }
            if e_cur[j] > v {
                v = e_cur[j];
                from = From::Left;
            }
            h_cur[j] = v;
            tb_h[i * w + j] = from;
            if v > best {
                (best, bi, bj) = (v, i, j);
            }
        }
        std::mem::swap(&mut h_prev, &mut h_cur);
        std::mem::swap(&mut f_prev, &mut f_cur);
    }
    if best <= 0 {
        return None;
    }
    // traceback
    let (mut i, mut j) = (bi, bj);
    let (mut pr, mut sr) = (Vec::new(), Vec::new());
    let mut state = tb_h[i * w + j];
    while i > 0 && j > 0 {
        match state {
            From::Stop => break,
            From::Diag => {
                pr.push(panel[i - 1].to_ascii_uppercase());
                sr.push(window[j - 1].to_ascii_uppercase());
                i -= 1;
                j -= 1;
                state = tb_h[i * w + j];
            }
            From::Up => {
                let ext = tb_f_ext[i * w + j];
                pr.push(panel[i - 1].to_ascii_uppercase());
                sr.push(b'-');
                i -= 1;
                state = if ext { From::Up } else { tb_h[i * w + j] };
            }
            From::Left => {
                let ext = tb_e_ext[i * w + j];
                pr.push(b'-');
                sr.push(window[j - 1].to_ascii_uppercase());
                j -= 1;
                state = if ext { From::Left } else { tb_h[i * w + j] };
            }
        }
    }
    pr.reverse();
    sr.reverse();
    let cols = pr.len() as u64;
    let same = pr
        .iter()
        .zip(&sr)
        .filter(|(a, b)| a == b && **a != b'-')
        .count() as u64;
    let gaps = pr
        .iter()
        .zip(&sr)
        .filter(|(a, b)| **a == b'-' || **b == b'-')
        .count() as u64;
    Some(PanelAlignment {
        panel_start: i as u64 + 1,
        panel_end: bi as u64,
        window_start: j as u64 + 1,
        window_end: bj as u64,
        identity: if cols > 0 {
            same as f64 * 100.0 / cols as f64
        } else {
            0.0
        },
        mismatches: cols - same - gaps,
        gap_bases: gaps,
        panel_row: pr,
        strain_row: sr,
    })
}

/// The strain stretch to align a panel gene against: the hit widened by
/// the gene's length each side (so a hit covering part of the gene can
/// still show the rest, or that it is missing), clipped to the contig,
/// in the hit's orientation. Returns (window, its first contig position).
pub fn hit_window(
    contig: &[u8],
    start: u64,
    end: u64,
    strand: i8,
    panel_len: usize,
) -> (Vec<u8>, u64) {
    let pad = panel_len.max(200) as u64;
    let lo = start.saturating_sub(pad).max(1);
    let hi = (end + pad).min(contig.len() as u64);
    let slice = &contig[(lo - 1) as usize..hi as usize];
    let w = if strand < 0 {
        revcomp(slice)
    } else {
        slice.to_ascii_uppercase()
    };
    (w, lo)
}

/// A window position (1-based, window orientation) on the contig.
pub fn window_to_contig(pos: u64, window_lo: u64, window_len: usize, strand: i8) -> u64 {
    if strand < 0 {
        window_lo + window_len as u64 - pos
    } else {
        window_lo + pos - 1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const GENE: &[u8] = b"ATGGCTAAAGAAACTGTTTATCGTGTAGACGGTTTATCTTGTACGAATTGTGCAGCTAAATTTGAG";

    #[test]
    fn finds_an_exact_copy_inside_flanks() {
        let window = [b"CCCCGGGGTTTT".as_slice(), GENE, b"AAAACCCC"].concat();
        let a = align(GENE, &window).unwrap();
        assert_eq!((a.panel_start, a.panel_end), (1, GENE.len() as u64));
        assert_eq!((a.window_start, a.window_end), (13, 12 + GENE.len() as u64));
        assert_eq!((a.identity, a.mismatches, a.gap_bases), (100.0, 0, 0));
        assert_eq!(a.panel_coverage(GENE.len()), 100.0);
    }

    #[test]
    fn shows_mismatches_and_a_deletion() {
        let mut copy = GENE.to_vec();
        copy[20] = b'A'; // T -> A
        copy.drain(40..43); // a codon lost
        let a = align(GENE, &copy).unwrap();
        assert_eq!(a.mismatches, 1);
        assert_eq!(a.gap_bases, 3);
        assert_eq!(a.strain_row.iter().filter(|&&c| c == b'-').count(), 3);
        assert_eq!(a.panel_row.len(), a.strain_row.len());
    }

    #[test]
    fn a_partial_copy_covers_part_of_the_gene() {
        let half = &GENE[..GENE.len() / 2];
        let a = align(GENE, half).unwrap();
        assert_eq!(a.panel_start, 1);
        assert!(a.panel_coverage(GENE.len()) < 60.0);
    }

    #[test]
    fn windows_follow_the_hit_strand() {
        let contig = b"AAAACCCCGGGGTTTT";
        let (w, lo) = hit_window(contig, 5, 8, -1, 2);
        assert_eq!(lo, 1);
        assert_eq!(w, revcomp(contig));
        assert_eq!(window_to_contig(1, lo, w.len(), -1), 16);
        assert_eq!(window_to_contig(1, lo, w.len(), 1), 1);
        assert!(align(b"ACGT", b"").is_none());
    }
}
