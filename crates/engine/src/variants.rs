//! Per-base variant events (SNPs, deletions, insertions) extracted from
//! the nucmer delta file, for the strain map markers and the
//! whole-genome alignment viewer.
//!
//! The events are sparse (strain vs strain: a few thousand positions
//! genome-wide), so they are cheap to ship to the viewer in full. They
//! are computed on demand from the kept delta files and cached next to
//! each query's result.json, so old runs work without re-running.

use crate::delta::{Alignment, DeltaFile};
use crate::fasta::{self, FastaRecord};
use crate::msa;
use crate::Result;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::path::Path;
use straincompass_types::{AlignmentEvents, DelEvent, InsEvent, SnpEvent};

/// Extract variant events per reference seqid for one query against the
/// reference, walking the reconstructed pairwise alignment of every
/// delta block.
///
/// Query bases are reported in reference orientation (the reconstruction
/// works on reverse complemented query slices for reverse-strand blocks),
/// so the viewer can compare them to the reference directly.
pub fn variant_events(
    delta: &DeltaFile,
    ref_records: &[FastaRecord],
    qry_records: &[FastaRecord],
) -> Result<BTreeMap<String, AlignmentEvents>> {
    let ref_by_id: HashMap<&str, &FastaRecord> =
        ref_records.iter().map(|r| (r.id.as_str(), r)).collect();
    let qry_by_id: HashMap<&str, &FastaRecord> =
        qry_records.iter().map(|r| (r.id.as_str(), r)).collect();

    // Per seqid accumulation. Overlapping nucmer blocks could report the
    // same event twice; the best (highest identity) block wins, so blocks
    // are processed in identity-descending order and claimed positions
    // are skipped afterwards.
    struct Acc {
        snps: Vec<SnpEvent>,
        dels: Vec<DelEvent>,
        ins: Vec<InsEvent>,
        claimed_snps: HashSet<u64>,
        claimed_dels: HashSet<u64>,
        claimed_ins: HashSet<u64>,
    }
    let mut per_seqid: BTreeMap<String, Acc> = BTreeMap::new();

    let mut order: Vec<&Alignment> = delta.alignments.iter().collect();
    order.sort_by(|a, b| {
        b.identity()
            .partial_cmp(&a.identity())
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    for a in order {
        let Some(ref_rec) = ref_by_id.get(a.ref_seqid.as_str()) else {
            continue;
        };
        let Some(qry_rec) = qry_by_id.get(a.qry_seqid.as_str()) else {
            continue;
        };
        let (ref_bases, qry_bases) = msa::block_bases(a, ref_rec, qry_rec);
        let pw = msa::reconstruct(a, &ref_bases, &qry_bases);
        let acc = per_seqid.entry(a.ref_seqid.clone()).or_insert_with(|| Acc {
            snps: Vec::new(),
            dels: Vec::new(),
            ins: Vec::new(),
            claimed_snps: HashSet::new(),
            claimed_dels: HashSet::new(),
            claimed_ins: HashSet::new(),
        });

        let mut c = 0;
        while c < pw.len() {
            match (pw.ref_row[c], pw.qry_row[c]) {
                (b'-', b'-') => c += 1,
                // insertion in the query: gap in the reference row.
                // ref_pos holds the position of the last emitted
                // reference base, i.e. the base the insertion follows.
                (b'-', _) => {
                    let pos = pw.ref_pos[c];
                    let mut seq: Vec<u8> = Vec::new();
                    while c < pw.len() && pw.ref_row[c] == b'-' && pw.qry_row[c] != b'-' {
                        seq.push(pw.qry_row[c]);
                        c += 1;
                    }
                    if !seq.is_empty() && acc.claimed_ins.insert(pos) {
                        acc.ins.push(InsEvent {
                            pos,
                            seq: String::from_utf8_lossy(&seq).into_owned(),
                        });
                    }
                }
                // insertion in the reference = deletion in the query.
                (_, b'-') => {
                    let pos = pw.ref_pos[c];
                    let mut len = 0u64;
                    while c < pw.len() && pw.qry_row[c] == b'-' && pw.ref_row[c] != b'-' {
                        len += 1;
                        c += 1;
                    }
                    if len > 0 && acc.claimed_dels.insert(pos) {
                        acc.dels.push(DelEvent { pos, len });
                    }
                }
                (r, q) => {
                    if r != q && acc.claimed_snps.insert(pw.ref_pos[c]) {
                        acc.snps.push(SnpEvent {
                            pos: pw.ref_pos[c],
                            r,
                            q,
                        });
                    }
                    c += 1;
                }
            }
        }
    }

    let mut out = BTreeMap::new();
    for (seqid, acc) in per_seqid {
        let mut events = AlignmentEvents {
            snps: acc.snps,
            dels: acc.dels,
            ins: acc.ins,
        };
        events.snps.sort_by_key(|e| e.pos);
        events.dels.sort_by_key(|e| e.pos);
        events.ins.sort_by_key(|e| e.pos);
        out.insert(seqid, events);
    }
    Ok(out)
}

/// File-based convenience: parse the FASTAs and the delta, extract the
/// events, write them as JSON to `dest` (atomic temp + rename).
pub fn write_variant_events(
    ref_fasta: &Path,
    qry_fasta: &Path,
    delta_path: &Path,
    dest: &Path,
) -> Result<()> {
    let ref_records = fasta::parse_fasta(ref_fasta)?;
    write_variant_events_with_ref(&ref_records, qry_fasta, delta_path, dest)
}

/// Same as [`write_variant_events`] but takes pre-parsed reference records,
/// so a caller computing events for several queries parses the reference
/// fasta (usually multi-megabase) only once.
pub fn write_variant_events_with_ref(
    ref_records: &[FastaRecord],
    qry_fasta: &Path,
    delta_path: &Path,
    dest: &Path,
) -> Result<()> {
    let qry_records = fasta::parse_fasta(qry_fasta)?;
    let delta = DeltaFile::parse(delta_path)?;
    let events = variant_events(&delta, ref_records, &qry_records)?;
    write_events_json(&events, dest)
}

fn write_events_json(events: &BTreeMap<String, AlignmentEvents>, dest: &Path) -> Result<()> {
    let bytes = serde_json::to_vec(events)
        .map_err(|e| crate::friendly(format!("Cannot encode the variant events: {e}")))?;
    if let Some(parent) = dest.parent() {
        std::fs::create_dir_all(parent)?;
    }
    let tmp = dest.with_extension("json.tmp");
    std::fs::write(&tmp, bytes)?;
    std::fs::rename(&tmp, dest)?;
    Ok(())
}
