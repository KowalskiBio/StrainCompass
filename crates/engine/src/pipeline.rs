//! Orchestration of one query-vs-reference comparison and its artifacts.

use crate::blast;
use crate::coverage;
use crate::delta::DeltaFile;
use crate::fasta;
use crate::gaps;
use crate::gff;
use crate::tools::ToolPaths;
use crate::Result;
use bactiment_types::{GapRow, GeneCoverageRow, PanelRow, RunParams, WgaBlock, WgaGene};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::io::Write;
use std::path::{Path, PathBuf};

/// Everything one comparison produces (one query against the reference).
#[derive(Debug, Serialize, Deserialize, Default)]
pub struct ComparisonResult {
    pub query_name: String,
    pub genes_coverage: Vec<GeneCoverageRow>,
    pub unaligned_gaps: Vec<GapRow>,
    pub panel: Option<Vec<PanelRow>>,
    /// Per seqid reference lengths from the reference fasta.
    pub ref_lengths: Vec<(String, u64)>,
    pub genes: Vec<WgaGene>,
    pub blocks: Vec<WgaBlock>,
    pub dnadiff_report: Option<String>,
}

/// Inputs needed for one comparison. Paths point to sanitized files
/// already staged by the caller.
pub struct ComparisonInputs<'a> {
    pub ref_fasta: &'a Path,
    pub ref_gff: &'a Path,
    pub qry_fasta: &'a Path,
    pub query_name: &'a str,
    pub panel_fasta: Option<&'a Path>,
    pub params: &'a RunParams,
}

/// Where the engine works and caches.
pub struct WorkDirs<'a> {
    /// scratch for this comparison
    pub work: &'a Path,
    /// persistent alignment cache (delta files survive across runs)
    pub cache: &'a Path,
}

/// Progress callback: (step message, done, total).
pub type Progress<'a> = &'a dyn Fn(&str, u32, u32);

/// Sanitize an uploaded fasta into `dest`, returning the sanitized path.
pub fn sanitize_into(source: &Path, dest: &Path) -> Result<()> {
    let recs = fasta::parse_fasta(source)?;
    fasta::write_sanitized_fasta(&recs, dest)?;
    Ok(())
}

/// Content hash of a file, for cache keys.
pub fn file_hash(path: &Path) -> std::io::Result<String> {
    let mut f = std::fs::File::open(path)?;
    let mut hasher = Sha256::new();
    let mut buf = [0u8; 65536];
    use std::io::Read;
    loop {
        let n = f.read(&mut buf)?;
        if n == 0 {
            break;
        }
        hasher.update(&buf[..n]);
    }
    Ok(format!("{:x}", hasher.finalize()))
}

/// The alignment cache key: reference content + query content + align
/// layer parameters.
pub fn align_cache_key(ref_hash: &str, qry_hash: &str, params: &RunParams) -> String {
    let mut h = Sha256::new();
    h.update(ref_hash.as_bytes());
    h.update(b"|");
    h.update(qry_hash.as_bytes());
    h.update(b"|");
    h.update(params.align_signature().as_bytes());
    format!("{:x}", h.finalize())
}

/// Run one full comparison. Delta alignment is reused from cache when
/// possible (cheap re-runs for postprocess parameter changes).
pub fn run_comparison(
    tools: &ToolPaths,
    inputs: &ComparisonInputs,
    dirs: &WorkDirs,
    progress: Progress,
) -> Result<ComparisonResult> {
    std::fs::create_dir_all(dirs.work)?;
    std::fs::create_dir_all(dirs.cache)?;

    progress("Reading the reference annotation", 0, 3);
    let genes = gff::parse_gff(inputs.ref_gff)?;
    let ref_records = fasta::parse_fasta(inputs.ref_fasta)?;
    let ref_lengths = fasta::seq_lengths(&ref_records);
    // Warn (silently in engine) if GFF seqids are missing from the fasta.
    let mut ref_lengths_vec: Vec<(String, u64)> = ref_lengths.clone().into_iter().collect();
    ref_lengths_vec.sort();

    let qry_records = fasta::parse_fasta(inputs.qry_fasta)?;

    progress("Aligning the query genome against the reference", 1, 3);
    let ref_hash = file_hash(inputs.ref_fasta)?;
    let qry_hash = file_hash(inputs.qry_fasta)?;
    let key = align_cache_key(&ref_hash, &qry_hash, inputs.params);
    let cache_delta = dirs.cache.join(format!("{key}.delta"));
    let (delta_path, cached) = tools.run_nucmer(
        inputs.ref_fasta,
        inputs.qry_fasta,
        dirs.work,
        "cmp",
        inputs.params.nucmer_minmatch,
        inputs.params.nucmer_breaklen,
        Some(&cache_delta),
    )?;
    let _ = cached;

    progress("Scoring genes by alignment coverage", 2, 3);
    let delta = DeltaFile::parse(&delta_path)?;
    let cov_rows = coverage::gene_coverage(
        &genes,
        &delta,
        &ref_records,
        &qry_records,
        inputs.params.present_cov,
        inputs.params.partial_cov,
    )?;
    let gap_rows = gaps::unaligned_gaps(&genes, &delta, &ref_lengths, inputs.params.min_gap);

    let panel = match inputs.panel_fasta {
        Some(p) => {
            let panel_dir = dirs.work.join("panel");
            let rows = blast::panel_recheck(
                tools,
                p,
                inputs.qry_fasta,
                &panel_dir,
                inputs.query_name,
                inputs.params.blast_cov,
                inputs.params.blast_pid,
                inputs.params.blast_evalue,
            )?;
            Some(rows)
        }
        None => None,
    };

    let dnadiff_report = if inputs.params.dnadiff {
        let cached_report = dirs.cache.join(format!("{key}.report"));
        if cached_report.is_file() {
            std::fs::read_to_string(&cached_report).ok()
        } else {
            let p = tools.run_dnadiff(inputs.ref_fasta, inputs.qry_fasta, dirs.work, "cmp")?;
            let text = std::fs::read_to_string(&p)?;
            let _ = std::fs::copy(&p, &cached_report);
            Some(text)
        }
    } else {
        None
    };

    let wga_genes: Vec<WgaGene> = genes
        .iter()
        .map(|g| WgaGene {
            locus_tag: g.locus_tag.clone(),
            symbol: g.symbol.clone(),
            biotype: g.biotype.clone(),
            seqid: g.seqid.clone(),
            start: g.start,
            end: g.end,
            strand: g.strand,
        })
        .collect();
    let blocks: Vec<WgaBlock> = delta
        .alignments
        .iter()
        .map(|a| WgaBlock {
            ref_seqid: a.ref_seqid.clone(),
            ref_start: a.ref_start,
            ref_end: a.ref_end,
            qry_seqid: a.qry_seqid.clone(),
            qry_start: a.qry_lo,
            qry_end: a.qry_hi,
            qry_rev: a.qry_rev,
            identity: a.identity(),
        })
        .collect();

    Ok(ComparisonResult {
        query_name: inputs.query_name.to_string(),
        genes_coverage: cov_rows,
        unaligned_gaps: gap_rows,
        panel,
        ref_lengths: ref_lengths_vec,
        genes: wga_genes,
        blocks,
        dnadiff_report,
    })
}

/// Write a genes coverage table as TSV.
pub fn write_genes_coverage_tsv(rows: &[GeneCoverageRow], out: &Path) -> std::io::Result<()> {
    let mut w = std::io::BufWriter::new(std::fs::File::create(out)?);
    writeln!(w, "locus_tag\tsymbol\tbiotype\tseqid\tstart\tend\tlength\tcov_bp\tcov_pct\tcall\tbest_identity\tmismatches\tindels")?;
    for r in rows {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{}\t{:.2}\t{}\t{}",
            r.locus_tag,
            r.symbol,
            r.biotype,
            r.seqid,
            r.start,
            r.end,
            r.length,
            r.cov_bp,
            r.cov_pct,
            r.call.as_str(),
            r.best_identity,
            r.mismatches,
            r.indels
        )?;
    }
    Ok(())
}

pub fn write_gaps_tsv(rows: &[GapRow], out: &Path) -> std::io::Result<()> {
    let mut w = std::io::BufWriter::new(std::fs::File::create(out)?);
    writeln!(w, "seqid\tstart\tend\tlength\tn_genes\tgenes")?;
    for r in rows {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}",
            r.seqid,
            r.start,
            r.end,
            r.length,
            r.n_genes,
            r.genes.join(";")
        )?;
    }
    Ok(())
}

pub fn write_panel_tsv(rows: &[PanelRow], out: &Path) -> std::io::Result<()> {
    let mut w = std::io::BufWriter::new(std::fs::File::create(out)?);
    writeln!(w, "gene_id\tqlen\tcov_pct\tidentity\tbest_evalue\tcall")?;
    for r in rows {
        writeln!(
            w,
            "{}\t{}\t{:.2}\t{:.2}\t{}\t{}",
            r.gene_id,
            r.qlen,
            r.cov_pct,
            r.identity,
            r.best_evalue,
            r.call.as_str()
        )?;
    }
    Ok(())
}

/// The run artifact paths, for the Files panel.
pub fn comparison_artifacts(run_dir: &Path, query_name: &str) -> Vec<(String, PathBuf)> {
    vec![
        (
            format!("{query_name}: genes coverage table"),
            run_dir.join("genes_coverage.tsv"),
        ),
        (
            format!("{query_name}: unaligned gaps table"),
            run_dir.join("unaligned_gaps.tsv"),
        ),
        (
            format!("{query_name}: alignment report (dnadiff)"),
            run_dir.join("dnadiff.report"),
        ),
    ]
}

/// Build the full gene detail (hover preview stats + MSA rows) for one
/// gene across all queries of a run. Computed on demand from the cached
/// delta artifacts; no re-alignment.
pub struct QueryAlignmentSource<'a> {
    pub query_id: i64,
    pub query_name: String,
    pub qry_fasta: &'a Path,
    pub delta: &'a Path,
}

pub fn gene_detail(
    ref_fasta: &Path,
    ref_gff: &Path,
    params: &RunParams,
    locus_tag: &str,
    sources: &[QueryAlignmentSource<'_>],
) -> Result<bactiment_types::GeneDetail> {
    let genes = gff::parse_gff(ref_gff)?;
    let gene = genes
        .iter()
        .find(|g| g.locus_tag == locus_tag)
        .ok_or_else(|| {
            crate::friendly(format!(
                "The gene {locus_tag} was not found in the reference annotation."
            ))
        })?;
    let ref_records = fasta::parse_fasta(ref_fasta)?;
    let ref_rec = ref_records
        .iter()
        .find(|r| r.id == gene.seqid)
        .ok_or_else(|| {
            crate::friendly(format!(
                "The reference sequence {} was not found.",
                gene.seqid
            ))
        })?;
    let gene_len = gene.end - gene.start + 1;
    let reference_seq = String::from_utf8_lossy(&fasta::subseq(
        ref_rec,
        gene.start,
        gene.end,
        gene.strand < 0,
    ))
    .into_owned();

    let mut queries = Vec::new();
    for src in sources {
        let qry_records = fasta::parse_fasta(src.qry_fasta)?;
        let delta = DeltaFile::parse(src.delta)?;
        // Coverage / call for this gene.
        let cov_rows = coverage::gene_coverage(
            std::slice::from_ref(gene),
            &delta,
            &ref_records,
            &qry_records,
            params.present_cov,
            params.partial_cov,
        )?;
        let cov = cov_rows.into_iter().next().unwrap_or_default();

        // Reconstruct blocks overlapping the gene.
        let mut blocks = Vec::new();
        for a in &delta.alignments {
            if a.ref_seqid != gene.seqid || a.ref_end < gene.start || a.ref_start > gene.end {
                continue;
            }
            let Some(qry_rec) = qry_records.iter().find(|r| r.id == a.qry_seqid) else {
                continue;
            };
            let (ref_bases, qry_bases) = crate::msa::block_bases(a, ref_rec, qry_rec);
            let pw = crate::msa::reconstruct(a, &ref_bases, &qry_bases);
            let slice = crate::msa::slice_block_to_ref_range(&pw, a, gene.start, gene.end);
            let (mut ref_row, mut qry_row) =
                crate::msa::orient_for_strand(gene.strand, slice.ref_row, slice.qry_row);
            if gene.strand < 0 {
                // qry positions were forward; display order stays fine
                // because we reversed the rows wholesale.
            }
            let _ = (&mut ref_row, &mut qry_row);
            blocks.push(bactiment_types::GeneBlock {
                ref_start: slice.ref_start.max(gene.start),
                ref_end: slice.ref_end.min(gene.end),
                qry_start: slice.qry_start,
                qry_end: slice.qry_end,
                qry_rev: a.qry_rev,
                identity: a.identity(),
                ref_seq: String::from_utf8_lossy(&ref_row).into_owned(),
                qry_seq: String::from_utf8_lossy(&qry_row).into_owned(),
            });
        }
        // Unaligned stretches inside the gene.
        let aligned = delta
            .aligned_ref_intervals()
            .get(gene.seqid.as_str())
            .cloned()
            .unwrap_or_default();
        let mut cov_intervals: Vec<(u64, u64)> = Vec::new();
        for (s, e) in aligned {
            let s = s.max(gene.start);
            let e = e.min(gene.end);
            if s <= e {
                cov_intervals.push((s, e));
            }
        }
        let mut unaligned = Vec::new();
        let mut cursor = gene.start;
        for (s, e) in cov_intervals {
            if s > cursor {
                unaligned.push((cursor, s - 1));
            }
            cursor = e + 1;
        }
        if cursor <= gene.end {
            unaligned.push((cursor, gene.end));
        }

        // Premature stops from the query sequence of this gene.
        let mut qry_seq_gene = Vec::new();
        for b in blocks.iter().flat_map(|b| b.qry_seq.bytes()) {
            if b != b'-' {
                qry_seq_gene.push(b);
            }
        }
        let stops = crate::msa::premature_stops(&qry_seq_gene)
            .into_iter()
            .map(|(ci, aa)| bactiment_types::PrematureStop {
                codon_index: ci,
                aa_position: aa,
            })
            .collect();

        queries.push(bactiment_types::GeneQueryAlignment {
            query_id: src.query_id,
            query_name: src.query_name.clone(),
            call: cov.call,
            cov_pct: cov.cov_pct,
            best_identity: cov.best_identity,
            mismatches: cov.mismatches,
            indels: cov.indels,
            blocks,
            unaligned,
            premature_stops: stops,
        });
    }

    Ok(bactiment_types::GeneDetail {
        locus_tag: gene.locus_tag.clone(),
        symbol: gene.symbol.clone(),
        biotype: gene.biotype.clone(),
        seqid: gene.seqid.clone(),
        start: gene.start,
        end: gene.end,
        strand: gene.strand,
        length: gene_len,
        reference_seq,
        queries,
    })
}
