//! Query-side unaligned regions ("gained" sequence), where they attach to
//! the reference, and the genes predicted inside them.
//!
//! This is the mirror image of [`crate::gaps`]: that module complements the
//! aligned intervals on the *reference* axis to find what a query lacks,
//! this one complements them on the *query* axis to find what a query has
//! and the reference does not.
//!
//! The asymmetry that matters: nucmer runs with its default
//! `--mumreference` matcher, whose anchors must be unique on the reference
//! side. A query copy of a gene the reference carries in several places
//! (an rRNA operon, an IS element, a transposase) may have no unique anchor
//! to seed from, fail to align, and land here even though the reference has
//! the gene three times over. So a region found here is one with *no
//! alignment to the reference*, which is weaker than *absent from the
//! reference*, and nothing downstream may claim otherwise.

use crate::delta::{Alignment, DeltaFile};
use crate::fasta::{self, FastaRecord};
use crate::gff::Gene;
use crate::tools::ToolPaths;
use crate::Result;
use std::collections::HashMap;
use std::path::Path;
use straincompass_types::{GainedAnchor, GainedOrfStatus, GainedRow};

/// How far apart the two flanks' reference positions may be and still count
/// as one insertion point. Beyond this the flanks are telling different
/// stories - a rearrangement, a repeat, a mis-assembly - and collapsing
/// them to a single position would be a lie.
const ANCHOR_MAX_SPAN: u64 = 20_000;

/// Compute the gained regions of one query: stretches of the query genome
/// (>= `min_gained` bp) that no alignment block covers, each placed on the
/// reference by the blocks flanking it, with the genes predicted inside it
/// when a gene finder is available.
///
/// Gene prediction is an enrichment on top of a result that is already
/// complete without it, so a missing or failing prodigal degrades to
/// [`GainedOrfStatus::Unavailable`] instead of failing the comparison. That
/// is a deliberate departure from blastn and dnadiff, which hard-fail:
/// those produce a result the user asked for, this only decorates one.
#[allow(clippy::too_many_arguments)]
pub fn gained_regions(
    tools: &ToolPaths,
    ref_fasta: &Path,
    ref_gff: &Path,
    delta: &DeltaFile,
    qry_records: &[FastaRecord],
    genes: &[Gene],
    region_fasta: &Path,
    work_dir: &Path,
    min_gained: u64,
    predict_orfs: bool,
) -> Result<(Vec<GainedRow>, GainedOrfStatus)> {
    let aligned = delta.aligned_qry_intervals();
    let qry_lengths = fasta::seq_lengths(qry_records);
    let by_id: HashMap<&str, &FastaRecord> =
        qry_records.iter().map(|r| (r.id.as_str(), r)).collect();

    // Raw blocks, not the merged intervals: merging has thrown away the
    // reference coordinates and the strand, which is exactly what anchoring
    // needs. Ascending by query start so the flank search can scan.
    let mut blocks_by_qry: HashMap<&str, Vec<&Alignment>> = HashMap::new();
    for a in &delta.alignments {
        blocks_by_qry
            .entry(a.qry_seqid.as_str())
            .or_default()
            .push(a);
    }
    for v in blocks_by_qry.values_mut() {
        v.sort_by_key(|a| (a.qry_lo, a.qry_hi));
    }

    // Reference genes sorted once per call, so each anchor resolves its
    // flanking genes by binary search.
    let mut genes_by_seqid: HashMap<&str, Vec<&Gene>> = HashMap::new();
    for g in genes {
        genes_by_seqid.entry(g.seqid.as_str()).or_default().push(g);
    }
    for v in genes_by_seqid.values_mut() {
        v.sort_by_key(|g| (g.start, g.end));
    }

    let mut rows: Vec<GainedRow> = Vec::new();
    let mut seqids: Vec<&String> = qry_lengths.keys().collect();
    seqids.sort();
    for seqid in seqids {
        let seqlen = qry_lengths[seqid];
        let ivs = aligned.get(seqid).cloned().unwrap_or_default();
        // Complement, exactly as gaps::unaligned_gaps does on the reference.
        let mut cursor = 1u64;
        let mut regions: Vec<(u64, u64)> = Vec::new();
        for (s, e) in &ivs {
            if *s > cursor {
                regions.push((cursor, s - 1));
            }
            cursor = e + 1;
        }
        if cursor <= seqlen {
            regions.push((cursor, seqlen));
        }

        let blocks = blocks_by_qry.get(seqid.as_str());
        for (start, end) in regions {
            let length = end - start + 1;
            if length < min_gained {
                continue;
            }
            let gc_pct = by_id
                .get(seqid.as_str())
                .map(|rec| gc_percent(&fasta::subseq(rec, start, end, false)))
                .unwrap_or(0.0);
            let mut row = GainedRow {
                qry_seqid: seqid.clone(),
                start,
                end,
                length,
                gc_pct,
                at_contig_end: start == 1 || end == seqlen,
                ..Default::default()
            };
            anchor_region(&mut row, blocks.map(|v| v.as_slice()), &genes_by_seqid);
            rows.push(row);
        }
    }
    rows.sort_by(|a, b| (a.qry_seqid.as_str(), a.start).cmp(&(b.qry_seqid.as_str(), b.start)));

    let status = predict_gained_orfs(
        tools,
        &mut rows,
        &by_id,
        region_fasta,
        work_dir,
        predict_orfs,
    )?;

    // Name what was predicted, in one translated search for the whole
    // query. Decoration on an already-complete result, like prediction
    // itself: a naming failure leaves the genes unnamed, it does not
    // take the regions down with it.
    if matches!(status, GainedOrfStatus::Predicted) {
        let _ =
            crate::blast::name_gained_orfs(tools, ref_fasta, ref_gff, &mut rows, &by_id, work_dir);
    }
    Ok((rows, status))
}

/// Place one region on the reference using the blocks that flank it in
/// query coordinates. Fills the anchor fields of `row` in place.
fn anchor_region(
    row: &mut GainedRow,
    blocks: Option<&[&Alignment]>,
    genes_by_seqid: &HashMap<&str, Vec<&Gene>>,
) {
    let Some(blocks) = blocks.filter(|b| !b.is_empty()) else {
        // No alignment anywhere on this contig: a whole extra replicon.
        row.anchor = GainedAnchor::Unanchored;
        return;
    };

    // Several blocks can end at the same merged boundary; prefer the longest,
    // which is the more trustworthy placement.
    let left = blocks
        .iter()
        .filter(|a| a.qry_hi < row.start)
        .max_by_key(|a| (a.qry_hi, a.qry_len()))
        .copied();
    let right = blocks
        .iter()
        .filter(|a| a.qry_lo > row.end)
        .min_by_key(|a| (a.qry_lo, std::cmp::Reverse(a.qry_len())))
        .copied();

    // On a reverse block the query's increasing coordinates run DOWN the
    // reference, so the query position just before the region corresponds to
    // the block's low reference coordinate, not its high one. Getting this
    // backwards puts the anchor off by the whole length of the block, and
    // does it silently.
    let left_ref = left.map(|a| if a.qry_rev { a.ref_start } else { a.ref_end });
    let right_ref = right.map(|a| if a.qry_rev { a.ref_end } else { a.ref_start });

    match (left, right) {
        (Some(l), Some(r)) => {
            let (lp, rp) = (left_ref.unwrap(), right_ref.unwrap());
            let far = lp.abs_diff(rp) > ANCHOR_MAX_SPAN;
            // Opposite orientations are not by themselves disqualifying: an
            // inversion boundary legitimately gives one forward and one
            // reverse flank still pointing at the same neighbourhood. The
            // distance is the test that matters; the rest is recorded.
            row.flanks_disagree = l.ref_seqid != r.ref_seqid || far || l.qry_rev != r.qry_rev;
            if l.ref_seqid == r.ref_seqid && !far {
                row.anchor = GainedAnchor::Between;
                row.anchor_seqid = l.ref_seqid.clone();
                row.anchor_start = lp.min(rp);
                row.anchor_end = lp.max(rp);
            } else {
                // Discordant: fall back to the left flank, deterministically.
                row.anchor = GainedAnchor::Flank;
                row.anchor_seqid = l.ref_seqid.clone();
                row.anchor_start = lp;
                row.anchor_end = lp;
            }
        }
        (Some(l), None) => {
            row.anchor = GainedAnchor::Flank;
            row.anchor_seqid = l.ref_seqid.clone();
            row.anchor_start = left_ref.unwrap();
            row.anchor_end = left_ref.unwrap();
        }
        (None, Some(r)) => {
            row.anchor = GainedAnchor::Flank;
            row.anchor_seqid = r.ref_seqid.clone();
            row.anchor_start = right_ref.unwrap();
            row.anchor_end = right_ref.unwrap();
        }
        (None, None) => {
            // Blocks exist on the contig, but none on either side of this
            // region - it spans everything that did align. Unplaceable.
            row.anchor = GainedAnchor::Unanchored;
            return;
        }
    }

    if let Some(gs) = genes_by_seqid.get(row.anchor_seqid.as_str()) {
        row.left_gene = gene_at_or_before(gs, row.anchor_start);
        row.right_gene = gene_at_or_after(gs, row.anchor_end);
    }
}

/// Locus tag of the last gene starting at or before `pos`: the gene
/// containing it, or the nearest one before it. Empty when there is none.
///
/// Genes are sorted by start, so this is one binary search - the region
/// count times a few thousand genes is not something to scan.
fn gene_at_or_before(genes: &[&Gene], pos: u64) -> String {
    let upper = genes.partition_point(|g| g.start <= pos);
    if upper == 0 {
        return String::new();
    }
    genes[upper - 1].locus_tag.clone()
}

/// Locus tag of the gene containing `pos`, else the nearest one starting
/// after it. Empty when there is none.
fn gene_at_or_after(genes: &[&Gene], pos: u64) -> String {
    let upper = genes.partition_point(|g| g.start <= pos);
    // The gene just before can still contain pos, in which case the
    // position is inside a gene rather than between two of them.
    if upper > 0 && genes[upper - 1].end >= pos {
        return genes[upper - 1].locus_tag.clone();
    }
    genes
        .get(upper)
        .map(|g| g.locus_tag.clone())
        .unwrap_or_default()
}

/// G+C over the unambiguous bases only. An N run is excluded from the
/// denominator rather than counted as not-GC, which would depress the
/// figure exactly where an assembly is least certain.
fn gc_percent(seq: &[u8]) -> f64 {
    let mut gc = 0u64;
    let mut acgt = 0u64;
    for b in seq {
        match b.to_ascii_uppercase() {
            b'G' | b'C' => {
                gc += 1;
                acgt += 1;
            }
            b'A' | b'T' => acgt += 1,
            _ => {}
        }
    }
    if acgt == 0 {
        0.0
    } else {
        100.0 * gc as f64 / acgt as f64
    }
}

/// Predict the genes inside each region with prodigal, filling `n_orfs`,
/// `n_orfs_complete` and `orfs` in place.
///
/// Returns the status rather than an error for every "could not predict"
/// case: the regions are already a complete result and must not be lost
/// because the server has no gene finder.
fn predict_gained_orfs(
    tools: &ToolPaths,
    rows: &mut [GainedRow],
    by_id: &HashMap<&str, &FastaRecord>,
    region_fasta: &Path,
    work_dir: &Path,
    predict: bool,
) -> Result<GainedOrfStatus> {
    if rows.is_empty() {
        // Nothing to look inside; saying so is more honest than "unavailable".
        return Ok(GainedOrfStatus::Predicted);
    }
    if !predict {
        return Ok(GainedOrfStatus::Unavailable(
            "Gene prediction was switched off for this run.".into(),
        ));
    }
    let Some(prodigal) = tools.prodigal.as_ref() else {
        return Ok(GainedOrfStatus::Unavailable(
            "The gene finder (prodigal) is not installed on this server, so the genes inside the gained regions could not be predicted. The regions themselves are unaffected.".into(),
        ));
    };

    // One record per region, named by its bare index. Anything richer in
    // the first token would have to be parsed back out, and every
    // separator we could pick can also occur in a sanitized contig id
    // (sanitize_id keeps '.', '_', '-' and ':'). The provenance goes in
    // the description instead, where prodigal ignores it and a human
    // downloading the file can still read it.
    write_region_fasta(rows, by_id, region_fasta)?;

    std::fs::create_dir_all(work_dir)?;
    let gff_out = work_dir.join("gained_orfs.gff");
    // -p meta: normal mode trains a GC and codon-usage model on its input
    // and balks below ~20 kb, which a gained set easily is. More to the
    // point, gained sequence is usually horizontally acquired - phage,
    // plasmid, IS - whose codon usage differs from the host, so training
    // on the host would be wrong and on the small biased gained set worse.
    // No -a/-d: protein and nucleotide output would roughly double this
    // feature's disk cost with nothing reading them yet. Adding a
    // "download predicted proteins" feature later is one flag away.
    let out = std::process::Command::new(prodigal)
        .args(["-p", "meta", "-f", "gff", "-q", "-i"])
        .arg(region_fasta)
        .arg("-o")
        .arg(&gff_out)
        .output();
    let out = match out {
        Ok(o) => o,
        Err(e) => {
            return Ok(GainedOrfStatus::Unavailable(format!(
                "The gene finder could not be started, so the genes inside the gained regions were not predicted. {e}"
            )))
        }
    };
    if !out.status.success() {
        return Ok(GainedOrfStatus::Unavailable(format!(
            "The gene finder did not finish correctly, so the genes inside the gained regions were not predicted. {}",
            String::from_utf8_lossy(&out.stderr).trim()
        )));
    }
    let text = match std::fs::read_to_string(&gff_out) {
        Ok(t) => t,
        Err(e) => {
            return Ok(GainedOrfStatus::Unavailable(format!(
                "The gene finder's output could not be read. {e}"
            )))
        }
    };

    // Every region gets a count now, including the ones with no genes:
    // from here on Some(0) means "we looked and found none".
    for r in rows.iter_mut() {
        r.orfs.clear();
        r.n_orfs = Some(0);
        r.n_orfs_complete = Some(0);
    }
    for (idx, orf) in parse_prodigal_gff(&text) {
        let Some(row) = rows.get_mut(idx) else {
            continue;
        };
        // Prodigal counts from 1 within the record it was given; shift to
        // query contig coordinates so an ORF can be pulled straight out of
        // query.fa without knowing which region it came from.
        let orf = straincompass_types::GainedOrf {
            start: row.start + orf.start - 1,
            end: row.start + orf.end - 1,
            ..orf
        };
        if !orf.partial {
            row.n_orfs_complete = Some(row.n_orfs_complete.unwrap_or(0) + 1);
        }
        row.n_orfs = Some(row.n_orfs.unwrap_or(0) + 1);
        row.orfs.push(orf);
    }
    Ok(GainedOrfStatus::Predicted)
}

/// Write one record per region, `>gr{index} {contig}:{start}-{end}`.
fn write_region_fasta(
    rows: &[GainedRow],
    by_id: &HashMap<&str, &FastaRecord>,
    out: &Path,
) -> Result<()> {
    use std::io::Write;
    if let Some(parent) = out.parent() {
        std::fs::create_dir_all(parent)?;
    }
    let mut w = std::io::BufWriter::new(std::fs::File::create(out)?);
    for (i, r) in rows.iter().enumerate() {
        let Some(rec) = by_id.get(r.qry_seqid.as_str()) else {
            continue;
        };
        let seq = fasta::subseq(rec, r.start, r.end, false);
        writeln!(w, ">gr{i} {}:{}-{}", r.qry_seqid, r.start, r.end)?;
        for chunk in seq.chunks(60) {
            w.write_all(chunk)?;
            w.write_all(b"\n")?;
        }
    }
    Ok(())
}

/// Parse prodigal's GFF into (region index, ORF) pairs, in region
/// coordinates. Unparseable lines are skipped rather than failing the
/// comparison - a gene finder's output is decoration, not a result.
fn parse_prodigal_gff(text: &str) -> Vec<(usize, straincompass_types::GainedOrf)> {
    let mut out = Vec::new();
    for line in text.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 9 {
            continue;
        }
        let Some(idx) = f[0]
            .strip_prefix("gr")
            .and_then(|v| v.parse::<usize>().ok())
        else {
            continue;
        };
        let (Ok(start), Ok(end)) = (f[3].parse::<u64>(), f[4].parse::<u64>()) else {
            continue;
        };
        let mut partial = false;
        let mut confidence = 0.0;
        for attr in f[8].split(';') {
            match attr.split_once('=') {
                // partial=00 means both ends are inside the region; anything
                // else means the ORF runs off an edge and is probably a
                // truncated gene rather than a whole one.
                Some(("partial", v)) => partial = v != "00",
                Some(("conf", v)) => confidence = v.parse().unwrap_or(0.0),
                _ => {}
            }
        }
        out.push((
            idx,
            straincompass_types::GainedOrf {
                start,
                end,
                strand: if f[6] == "-" { -1 } else { 1 },
                partial,
                confidence,
                best: None,
            },
        ));
    }
    out
}
