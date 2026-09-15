//! End-to-end engine test against real MUMmer + BLAST+.
//! Skipped when the tools are not on this machine (e.g. plain CI).

use bactiment_engine::pipeline::{self, ComparisonInputs, WorkDirs};
use bactiment_engine::tools::ToolPaths;
use bactiment_types::{Call, RunParams};
use std::io::Write;

fn tools_available() -> Option<ToolPaths> {
    ToolPaths::discover().ok()
}

struct XorShift(u64);
impl XorShift {
    fn next(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x
    }
    fn base(&mut self) -> u8 {
        b"ACGT"[(self.next() % 4) as usize]
    }
}

fn write_fasta(path: &std::path::Path, id: &str, seq: &[u8]) {
    let mut f = std::fs::File::create(path).unwrap();
    writeln!(f, ">{id} synthetic").unwrap();
    for chunk in seq.chunks(60) {
        writeln!(f, "{}", String::from_utf8_lossy(chunk)).unwrap();
    }
}

/// Genes G0..G7: 990bp each, starts at 2000 + i*2000 (1-based),
/// alternating strand (even = plus), G3 is tRNA.
fn genes() -> Vec<(u64, u64, i8, &'static str, &'static str)> {
    (0..8u64)
        .map(|i| {
            let start = 2000 + i * 2000;
            let end = start + 989;
            let strand = if i % 2 == 0 { 1i8 } else { -1 };
            let tag: &'static str = Box::leak(format!("G{i}").into_boxed_str());
            let biotype = if i == 3 { "tRNA" } else { "protein_coding" };
            (start, end, strand, tag, biotype)
        })
        .collect()
}

fn gene_start(i: u64) -> u64 {
    2000 + i * 2000
}
fn gene_end(i: u64) -> u64 {
    gene_start(i) + 989
}

fn write_gff(path: &std::path::Path) {
    let mut f = std::fs::File::create(path).unwrap();
    writeln!(f, "##gff-version 3").unwrap();
    for (i, (start, end, strand, tag, biotype)) in genes().iter().enumerate() {
        writeln!(
            f,
            "chr1\tsyn\tgene\t{start}\t{end}\t.\t{}\t.\tID=gene_{i};locus_tag={tag};gene=sym_{i};gene_biotype={biotype}",
            if *strand > 0 { "+" } else { "-" }
        )
        .unwrap();
    }
}

/// Reference: 40 kb random. Query derived from it:
/// - SNP inside G0 (plus strand) at gene offset 3
/// - 10bp insertion inside G1 at gene offset 500
/// - 1489bp deletion removing G3 entirely
/// - premature stop (TTA on the reference => TAA on the minus strand of
///   G5) at a codon boundary of G5
/// - query truncated inside G6 so G6 is PARTIAL
fn build_fixture(dir: &std::path::Path) -> Vec<u8> {
    let mut rng = XorShift(0xBADC0DE);
    let mut ref_seq: Vec<u8> = (0..40_000).map(|_| rng.base()).collect();
    // ATG at the start of G0 (plus strand) for realistic translation.
    let g0 = gene_start(0);
    ref_seq[(g0 - 1) as usize..(g0 + 2) as usize].copy_from_slice(b"ATG");
    // Premature stop inside G5 (minus strand): the gene strand reads the
    // revcomp, so a TTA edit on the reference strand becomes TAA (a stop).
    // Pin the reference bases to CCC so the edit differs at all three
    // positions. The codon at ref 0-based 12300..12303 is a gene codon
    // boundary (gene_end(5) - 12302 = 687 = 3 * 229).
    let stop_pos = (gene_start(5) + 300) as usize;
    ref_seq[stop_pos..stop_pos + 3].copy_from_slice(b"CCC");

    write_gff(&dir.join("ref.gff"));
    write_fasta(&dir.join("ref.fa"), "chr1", &ref_seq);

    // ---- build the query ----
    let mut q: Vec<u8> = ref_seq.clone();
    q[stop_pos..stop_pos + 3].copy_from_slice(b"TTA");
    // SNP inside G0 (ensure it differs too)
    {
        let pos = g0 as usize + 3;
        let orig = q[pos];
        q[pos] = match orig {
            b'A' => b'C',
            b'C' => b'G',
            b'G' => b'T',
            _ => b'A',
        };
    }

    // 10bp insertion inside G1
    let ins_at = (gene_start(1) + 500) as usize;
    let mut q2: Vec<u8> = q[..ins_at].to_vec();
    q2.extend(std::iter::repeat_n(b'A', 10));
    q2.extend_from_slice(&q[ins_at..]);

    // 1489bp deletion removing G3 entirely (plus flanks)
    let del_from = (gene_start(3) - 500) as usize;
    let del_to = (gene_end(3) + 500) as usize;
    let mut q3: Vec<u8> = q2[..del_from].to_vec();
    q3.extend_from_slice(&q2[del_to..]);

    // Truncate inside G6: cover roughly half of the gene.
    // The deleted span is q2[del_from..del_to] = 1989 bases and the
    // insertion added 10, so for ref 0-based r >= del_to the q3 index is
    // r - 1979. Cut right after ref 0-based 14499 (= 1-based 14500).
    let cut_ref = (gene_start(6) + 500 - 1) as usize; // ref 0-based 14499
    let cut_q = cut_ref - 1979 + 1;
    let q_final: Vec<u8> = q3[..cut_q].to_vec();
    write_fasta(&dir.join("query.fa"), "ctg1", &q_final);
    ref_seq
}

fn params() -> RunParams {
    RunParams {
        min_gap: 200,
        present_cov: 95.0,
        partial_cov: 1.0,
        ..RunParams::default()
    }
}

#[test]
fn full_pipeline_with_real_tools() {
    let Some(tools) = tools_available() else {
        eprintln!("skipping: MUMmer/BLAST+ not found");
        return;
    };
    let dir = std::env::temp_dir().join(format!("bactiment-test-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let ref_seq = build_fixture(&dir);

    let params = params();
    // Panel: reference slices of G0 and G3.
    {
        let g0 = &ref_seq[(gene_start(0) - 1) as usize..gene_end(0) as usize];
        let g3 = &ref_seq[(gene_start(3) - 1) as usize..gene_end(3) as usize];
        let mut f = std::fs::File::create(dir.join("panel.fa")).unwrap();
        writeln!(f, ">G0").unwrap();
        writeln!(f, "{}", String::from_utf8_lossy(g0)).unwrap();
        writeln!(f, ">G3").unwrap();
        writeln!(f, "{}", String::from_utf8_lossy(g3)).unwrap();
    }

    let inputs = ComparisonInputs {
        ref_fasta: &dir.join("ref.fa"),
        ref_gff: &dir.join("ref.gff"),
        qry_fasta: &dir.join("query.fa"),
        query_name: "testq",
        panel_fasta: Some(&dir.join("panel.fa")),
        params: &params,
    };
    let work = dir.join("work");
    let dirs = WorkDirs {
        work: &work,
        cache: &dir.join("cache"),
    };
    let progress = |_m: &str, _d: u32, _t: u32| {};
    let res = pipeline::run_comparison(&tools, &inputs, &dirs, &progress).unwrap();

    let by_tag = |t: &str| {
        res.genes_coverage
            .iter()
            .find(|r| r.locus_tag == t)
            .unwrap_or_else(|| panic!("gene {t} missing"))
            .clone()
    };

    let g0 = by_tag("G0");
    assert_eq!(g0.call, Call::Present);
    assert_eq!(g0.mismatches, 1, "SNP inside G0, got {:?}", g0);
    let g1 = by_tag("G1");
    assert_eq!(g1.call, Call::Present);
    assert_eq!(g1.indels, 10, "10bp insertion inside G1");
    let g3 = by_tag("G3");
    assert_eq!(g3.call, Call::Absent, "G3 removed by the deletion");
    assert_eq!(g3.cov_bp, 0);
    let g6 = by_tag("G6");
    assert_eq!(
        g6.call,
        Call::Partial,
        "G6 truncated at contig end, got {:?}",
        g6
    );

    // Gaps: the deletion appears as a ~1.99kb unaligned region containing
    // G3, and the truncation leaves the rest of the reference unaligned.
    let del_gap = res
        .unaligned_gaps
        .iter()
        .find(|g| (1900..2100).contains(&g.length))
        .expect("the deletion gap");
    assert!(
        del_gap.genes.contains(&"G3".to_string()),
        "G3 inside the gap"
    );
    let _end_gap = res
        .unaligned_gaps
        .iter()
        .find(|g| g.start == 14501 && g.length > 25000)
        .expect("the truncation gap");

    // Panel recheck: G0 present, G3 absent
    let panel = res.panel.as_ref().unwrap();
    let pg0 = panel.iter().find(|r| r.gene_id == "G0").unwrap();
    let pg3 = panel.iter().find(|r| r.gene_id == "G3").unwrap();
    assert_eq!(pg0.call, Call::Present);
    assert_eq!(pg3.call, Call::Absent);

    // dnadiff report produced
    let report = res.dnadiff_report.as_ref().unwrap();
    assert!(report.contains("TotalBases"));

    // ---- gene detail (MSA viewer data) for the minus strand gene G5 ----
    let qry_fa = dir.join("query.fa");
    let delta_path = work.join("cmp.delta");
    let ref_fa = dir.join("ref.fa");
    let ref_gff = dir.join("ref.gff");
    let sources = vec![pipeline::QueryAlignmentSource {
        query_id: 1,
        query_name: "testq".into(),
        qry_fasta: &qry_fa,
        delta: &delta_path,
    }];
    let detail = pipeline::gene_detail(&ref_fa, &ref_gff, &params, "G5", &sources).unwrap();
    assert_eq!(detail.queries.len(), 1);
    let q = &detail.queries[0];
    assert_eq!(q.call, Call::Present);
    assert!(
        !q.premature_stops.is_empty(),
        "expected a premature stop in G5"
    );
    assert_eq!(detail.reference_seq.len(), 990);
    let joined_ref: Vec<u8> = q.blocks.iter().flat_map(|b| b.ref_seq.bytes()).collect();
    assert_eq!(joined_ref.len(), 990, "G5 fully covered by blocks");
    let joined_qry: Vec<u8> = q
        .blocks
        .iter()
        .flat_map(|b| b.qry_seq.bytes())
        .filter(|&c| c != b'-')
        .collect();
    let diff = joined_qry
        .iter()
        .zip(joined_ref.iter())
        .filter(|(a, b)| a != b)
        .count();
    assert_eq!(
        diff, 3,
        "the TTA edit appears as 3 mismatches after revcomp, got {diff}"
    );

    // ---- cache: second run with only postprocess changes ----
    let params2 = RunParams {
        present_cov: 50.0,
        ..params.clone()
    };
    let inputs2 = ComparisonInputs {
        ref_fasta: &dir.join("ref.fa"),
        ref_gff: &dir.join("ref.gff"),
        qry_fasta: &dir.join("query.fa"),
        query_name: "testq",
        panel_fasta: None,
        params: &params2,
    };
    let res2 = pipeline::run_comparison(&tools, &inputs2, &dirs, &progress).unwrap();
    let g6b = res2
        .genes_coverage
        .iter()
        .find(|r| r.locus_tag == "G6")
        .unwrap();
    assert_eq!(g6b.call, Call::Present);

    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn reverse_complement_query_full_coverage() {
    let Some(tools) = tools_available() else {
        eprintln!("skipping: MUMmer/BLAST+ not found");
        return;
    };
    let dir = std::env::temp_dir().join(format!("bactiment-rc-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let ref_seq = build_fixture(&dir);
    let rc = bactiment_engine::fasta::revcomp(&ref_seq);
    write_fasta(&dir.join("qrc.fa"), "rc_ctg", &rc);

    let params = params();
    let inputs = ComparisonInputs {
        ref_fasta: &dir.join("ref.fa"),
        ref_gff: &dir.join("ref.gff"),
        qry_fasta: &dir.join("qrc.fa"),
        query_name: "rc",
        panel_fasta: None,
        params: &params,
    };
    let work = dir.join("work");
    let dirs = WorkDirs {
        work: &work,
        cache: &dir.join("cache"),
    };
    let progress = |_m: &str, _d: u32, _t: u32| {};
    let res = pipeline::run_comparison(&tools, &inputs, &dirs, &progress).unwrap();
    for row in &res.genes_coverage {
        assert_eq!(
            row.call,
            Call::Present,
            "gene {} not present",
            row.locus_tag
        );
        assert!(row.cov_pct > 99.99);
    }
    assert!(res.unaligned_gaps.iter().all(|g| g.length < 200));
    assert!(res.blocks.iter().all(|b| b.qry_rev));

    // MSA of a minus strand gene must reconstruct exactly
    let qry_fa = dir.join("qrc.fa");
    let delta_path = work.join("cmp.delta");
    let ref_fa = dir.join("ref.fa");
    let ref_gff = dir.join("ref.gff");
    let sources = vec![pipeline::QueryAlignmentSource {
        query_id: 1,
        query_name: "rc".into(),
        qry_fasta: &qry_fa,
        delta: &delta_path,
    }];
    let detail = pipeline::gene_detail(&ref_fa, &ref_gff, &params, "G0", &sources).unwrap();
    assert_eq!(detail.reference_seq.len(), 990);
    let q = &detail.queries[0];
    let joined_ref: Vec<u8> = q.blocks.iter().flat_map(|b| b.ref_seq.bytes()).collect();
    assert_eq!(joined_ref.len(), 990);
    let mismatches: usize = q
        .blocks
        .iter()
        .flat_map(|b| b.ref_seq.bytes().zip(b.qry_seq.bytes()))
        .filter(|(r, s)| *r != b'-' && *s != b'-' && r != s)
        .count();
    assert_eq!(
        mismatches, 0,
        "revcomp of the reference must match perfectly"
    );
    let _ = std::fs::remove_dir_all(&dir);
}
