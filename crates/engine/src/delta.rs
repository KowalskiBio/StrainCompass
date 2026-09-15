//! MUMmer delta file parsing and alignment reconstruction.
//!
//! The delta format (nucmer output), verified empirically against
//! show-coords / show-aligns / show-snps:
//!
//! ```text
//! >reference_file query_file
//! reflen qrylen
//! s1 e1 s2 e2 errors errors2 sim
//! <signed delta values...>
//! 0
//! ```
//!
//! Conventions (locked by unit tests against real nucmer output):
//! - `s1 e1` are ascending reference coordinates, `s2 e2` query
//!   coordinates; `s2 > e2` means the query aligns reverse-complemented.
//! - `errors` (5th field) = mismatch count + total indel bases.
//! - Each delta value encodes ONE gap column: advance both sequences by
//!   `|d| - 1` paired bases, then a single gap column follows on the side
//!   given by the sign: `d < 0` = insertion in the query (gap in the
//!   reference row), `d > 0` = insertion in the reference (gap in the
//!   query row).
//! - Identity (as show-coords reports it) = matches / total columns,
//!   where total columns = paired bases + all gap columns.

use crate::{friendly, Result};
use std::collections::HashMap;
use std::io::Read;
use std::path::Path;

#[derive(Debug, Clone)]
pub struct Alignment {
    pub ref_seqid: String,
    pub qry_seqid: String,
    /// Ascending reference coordinates (1-based inclusive).
    pub ref_start: u64,
    pub ref_end: u64,
    /// Query coordinates (1-based inclusive, forward orientation).
    pub qry_lo: u64,
    pub qry_hi: u64,
    pub qry_rev: bool,
    /// 5th header field: mismatches + indel bases.
    pub errors: u64,
    /// Signed delta values; each is one gap column.
    pub deltas: Vec<i64>,
}

impl Alignment {
    pub fn ref_len(&self) -> u64 {
        self.ref_end - self.ref_start + 1
    }
    pub fn qry_len(&self) -> u64 {
        self.qry_hi - self.qry_lo + 1
    }
    /// Bases present in the query but not in the reference (insertions in
    /// the query; they appear as dots in the reference row). Every delta
    /// value encodes exactly one gap column.
    pub fn qry_ins(&self) -> u64 {
        self.deltas.iter().filter(|d| **d < 0).count() as u64
    }
    /// Bases present in the reference but not in the query (insertions in
    /// the reference; dots in the query row).
    pub fn ref_ins(&self) -> u64 {
        self.deltas.iter().filter(|d| **d > 0).count() as u64
    }
    pub fn indels(&self) -> u64 {
        self.qry_ins() + self.ref_ins()
    }
    /// Columns where both sequences have a base.
    pub fn paired(&self) -> u64 {
        self.ref_len().saturating_sub(self.ref_ins())
    }
    pub fn mismatches(&self) -> u64 {
        self.errors.saturating_sub(self.indels())
    }
    pub fn matches(&self) -> u64 {
        self.paired().saturating_sub(self.mismatches())
    }
    /// Total alignment columns (what show-aligns prints).
    pub fn aligned_len(&self) -> u64 {
        self.paired() + self.indels()
    }
    /// Percent identity over all alignment columns (show-coords semantics).
    pub fn identity(&self) -> f64 {
        let cols = self.aligned_len();
        if cols == 0 {
            return 0.0;
        }
        100.0 * self.matches() as f64 / cols as f64
    }
}

#[derive(Debug, Clone, Default)]
pub struct DeltaFile {
    pub alignments: Vec<Alignment>,
}

impl DeltaFile {
    pub fn parse<P: AsRef<Path>>(path: P) -> Result<DeltaFile> {
        let mut text = String::new();
        std::fs::File::open(path.as_ref())?.read_to_string(&mut text)?;
        parse_delta_str(&text)
    }

    /// Merged aligned intervals on the reference: seqid -> sorted,
    /// non-overlapping intervals.
    pub fn aligned_ref_intervals(&self) -> HashMap<String, Vec<(u64, u64)>> {
        let mut map: HashMap<String, Vec<(u64, u64)>> = HashMap::new();
        for a in &self.alignments {
            map.entry(a.ref_seqid.clone())
                .or_default()
                .push((a.ref_start, a.ref_end));
        }
        for ivs in map.values_mut() {
            merge_intervals(ivs);
        }
        map
    }
}

/// Sort and merge overlapping/adjacent intervals in place.
pub fn merge_intervals(ivs: &mut Vec<(u64, u64)>) {
    if ivs.is_empty() {
        return;
    }
    ivs.sort();
    let mut merged: Vec<(u64, u64)> = Vec::with_capacity(ivs.len());
    merged.push(ivs[0]);
    for &(s, e) in ivs.iter().skip(1) {
        let last = merged.last_mut().unwrap();
        if s <= last.1.saturating_add(1) {
            last.1 = last.1.max(e);
        } else {
            merged.push((s, e));
        }
    }
    *ivs = merged;
}

pub fn parse_delta_str(text: &str) -> Result<DeltaFile> {
    let mut alignments = Vec::new();
    let mut lines = text.lines().peekable();
    // Skip the file-level header lines (paths + "NUCMER") until the first
    // sequence-pair header.
    while let Some(line) = lines.peek() {
        if line.trim().starts_with('>') {
            break;
        }
        lines.next();
    }

    let mut current: Option<(String, String)> = None;
    while let Some(line) = lines.next() {
        let line = line.trim_end_matches(['\r', '\n']);
        if line.is_empty() {
            continue;
        }
        if let Some(header) = line.strip_prefix('>') {
            // ">refid qryid reflen qrylen"
            let mut parts = header.split_whitespace();
            let ref_seqid = parts.next().unwrap_or("").to_string();
            let qry_seqid = parts.next().unwrap_or("").to_string();
            current = Some((ref_seqid, qry_seqid));
            continue;
        }
        let (ref_seqid, qry_seqid) = match &current {
            Some(c) => c.clone(),
            None => {
                return Err(friendly(
                    "The alignment file (delta) is malformed: data before header.",
                ));
            }
        };
        let fields: Vec<i64> = line
            .split_whitespace()
            .map(|f| f.parse::<i64>())
            .collect::<std::result::Result<_, _>>()
            .map_err(|_| {
                friendly("The alignment file (delta) is malformed: bad alignment header.")
            })?;
        if fields.len() != 7 {
            return Err(friendly(
                "The alignment file (delta) is malformed: expected 7 numbers in the alignment header.",
            ));
        }
        let (s1, e1, s2, e2, errors) = (fields[0], fields[1], fields[2], fields[3], fields[4]);
        let mut deltas = Vec::new();
        for dline in lines.by_ref() {
            let dline = dline.trim();
            if dline.is_empty() {
                continue;
            }
            match dline.parse::<i64>() {
                Ok(0) => break,
                Ok(v) => deltas.push(v),
                Err(_) => {
                    return Err(friendly(
                        "The alignment file (delta) is malformed: bad delta value.",
                    ));
                }
            }
        }
        let qry_rev = e2 < s2;
        let qry_lo = s2.min(e2).unsigned_abs();
        let qry_hi = s2.max(e2).unsigned_abs();
        alignments.push(Alignment {
            ref_seqid,
            qry_seqid,
            ref_start: s1.unsigned_abs(),
            ref_end: e1.unsigned_abs(),
            qry_lo,
            qry_hi,
            qry_rev,
            errors: errors.unsigned_abs(),
            deltas,
        });
    }
    Ok(DeltaFile { alignments })
}
