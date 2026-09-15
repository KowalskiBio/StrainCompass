//! FASTA reading, writing, sanitization and subsequence extraction.

use crate::{friendly, Result};
use std::collections::HashMap;
use std::io::{BufWriter, Read, Write};
use std::path::Path;

#[derive(Debug, Clone)]
pub struct FastaRecord {
    pub id: String,
    pub desc: String,
    pub seq: Vec<u8>,
}

impl FastaRecord {
    pub fn seq_str(&self) -> String {
        String::from_utf8_lossy(&self.seq).into_owned()
    }
}

/// Parse a FASTA file. Validatates that the content looks like FASTA and
/// that sequence characters are sane, with plain language errors.
pub fn parse_fasta<P: AsRef<Path>>(path: P) -> Result<Vec<FastaRecord>> {
    let mut text = String::new();
    std::fs::File::open(path.as_ref())?.read_to_string(&mut text)?;
    parse_fasta_str(&text)
}

pub fn parse_fasta_str(text: &str) -> Result<Vec<FastaRecord>> {
    let mut records = Vec::new();
    let mut cur_id: Option<String> = None;
    let mut cur_desc = String::new();
    let mut cur_seq: Vec<u8> = Vec::new();
    let mut seen_any_line = false;
    for line in text.lines() {
        let line = line.trim_end_matches(['\r', '\n']);
        if line.is_empty() {
            continue;
        }
        seen_any_line = true;
        if let Some(header) = line.strip_prefix('>') {
            if let Some(id) = cur_id.take() {
                records.push(FastaRecord {
                    id,
                    desc: std::mem::take(&mut cur_desc),
                    seq: std::mem::take(&mut cur_seq),
                });
            }
            let header = header.trim();
            if header.is_empty() {
                return Err(friendly(
                    "This FASTA file has a sequence with an empty name. Please check the file.",
                ));
            }
            let mut parts = header.split_whitespace();
            cur_id = Some(parts.next().unwrap().to_string());
            cur_desc = parts.collect::<Vec<_>>().join(" ");
        } else {
            if cur_id.is_none() {
                return Err(friendly(
                    "This does not look like a FASTA file. A FASTA file starts with a line like \">sequence name\".",
                ));
            }
            for b in line.bytes() {
                let b = b.to_ascii_uppercase();
                match b {
                    b'A' | b'C' | b'G' | b'T' | b'U' | b'R' | b'Y' | b'S' | b'W' | b'K' | b'M'
                    | b'B' | b'D' | b'H' | b'V' | b'N' | b'-' | b'*' => cur_seq.push(b),
                    _ => {
                        return Err(friendly(format!(
                            "This FASTA file contains an unexpected character ({}) in the sequence. The file may be corrupted or not a FASTA file.",
                            b as char
                        )));
                    }
                }
            }
        }
    }
    if !seen_any_line {
        return Err(friendly("This file is empty."));
    }
    if let Some(id) = cur_id {
        records.push(FastaRecord {
            id,
            desc: cur_desc,
            seq: cur_seq,
        });
    }
    if records.is_empty() {
        return Err(friendly("No sequences were found in this FASTA file."));
    }
    Ok(records)
}

/// Replace a sequence id with a safe id: no whitespace, only
/// [A-Za-z0-9._:-]. Keeps a mapping so outputs can use original names.
pub fn sanitize_id(raw: &str) -> String {
    let mut out = String::with_capacity(raw.len());
    for c in raw.chars() {
        if c.is_ascii_alphanumeric() || matches!(c, '.' | '_' | '-' | ':') {
            out.push(c);
        } else {
            out.push('_');
        }
    }
    if out.is_empty() {
        out.push('s');
    }
    // Avoid leading dashes which some tools dislike.
    if out.starts_with('-') {
        out.replace_range(0..1, "_");
    }
    out
}

/// Write a sanitized copy of the input FASTA: safe unique ids, uppercase
/// sequences, no description. Returns (records, map sanitized -> original).
pub fn write_sanitized_fasta<P: AsRef<Path>>(
    source: &[FastaRecord],
    dest: P,
) -> Result<HashMap<String, String>> {
    let mut map = HashMap::new();
    let mut used: HashMap<String, usize> = HashMap::new();
    let file = std::fs::File::create(dest.as_ref())?;
    let mut w = BufWriter::new(file);
    for rec in source {
        let mut id = sanitize_id(&rec.id);
        let n = used.entry(id.clone()).or_insert(0);
        if *n > 0 {
            id = format!("{id}_{n}");
        }
        *n += 1;
        map.insert(id.clone(), rec.id.clone());
        writeln!(w, ">{id}")?;
        let seq = rec.seq.to_ascii_uppercase();
        for chunk in seq.chunks(60) {
            writeln!(w, "{}", String::from_utf8_lossy(chunk))?;
        }
    }
    w.flush()?;
    Ok(map)
}

pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|b| match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'U' => b'A',
            other => complement_iupac(other),
        })
        .collect()
}

fn complement_iupac(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'R' => b'Y',
        b'Y' => b'R',
        b'S' => b'S',
        b'W' => b'W',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        _ => b'N',
    }
}

/// Extract a 1-based inclusive range from a record; handles strand.
pub fn subseq(rec: &FastaRecord, start: u64, end: u64, rev: bool) -> Vec<u8> {
    let s = (start.max(1) as usize).saturating_sub(1);
    let e = (end as usize).min(rec.seq.len());
    let mut slice = rec.seq[s..e.max(s)].to_vec();
    if rev {
        slice = revcomp(&slice);
    }
    slice
}

/// Contig lengths by id.
pub fn seq_lengths(records: &[FastaRecord]) -> HashMap<String, u64> {
    records
        .iter()
        .map(|r| (r.id.clone(), r.seq.len() as u64))
        .collect()
}
