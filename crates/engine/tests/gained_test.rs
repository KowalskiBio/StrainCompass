use std::path::PathBuf;
use straincompass_engine::delta::parse_delta_str;
use straincompass_engine::fasta::FastaRecord;
use straincompass_engine::gained::gained_regions;
use straincompass_engine::gff::Gene;
use straincompass_engine::tools::ToolPaths;
use straincompass_types::{GainedAnchor, GainedOrfStatus};

/// A delta header line plus its terminating 0, for one gap-free block.
/// `s2`/`e2` are written as given, so passing e2 < s2 makes it reverse.
fn block(
    ref_id: &str,
    qry_id: &str,
    reflen: u64,
    qrylen: u64,
    coords: &[(i64, i64, i64, i64)],
) -> String {
    let mut s = format!("/ref.fa /qry.fa\nNUCMER\n>{ref_id} {qry_id} {reflen} {qrylen}\n");
    for (s1, e1, s2, e2) in coords {
        s.push_str(&format!("{s1} {e1} {s2} {e2} 0 0 0\n0\n"));
    }
    s
}

fn rec(id: &str, seq: &str) -> FastaRecord {
    FastaRecord {
        id: id.into(),
        desc: String::new(),
        seq: seq.as_bytes().to_vec(),
    }
}

/// A query contig of `len` bp, all A so GC is 0 unless a test says otherwise.
fn flat(id: &str, len: usize) -> FastaRecord {
    rec(id, &"A".repeat(len))
}

fn gene(seqid: &str, locus: &str, start: u64, end: u64) -> Gene {
    Gene {
        seqid: seqid.into(),
        start,
        end,
        strand: 1,
        locus_tag: locus.into(),
        old_locus_tag: String::new(),
        symbol: String::new(),
        biotype: "protein_coding".into(),
        protein_id: String::new(),
        product: String::new(),
    }
}

/// Tools with no gene finder: Phase 1 never predicts, and every real
/// deployment that lacks prodigal takes this path too.
fn no_tools() -> ToolPaths {
    let p = PathBuf::from("/nonexistent");
    ToolPaths {
        nucmer: p.clone(),
        show_coords: p.clone(),
        show_snps: p.clone(),
        dnadiff: p.clone(),
        makeblastdb: p.clone(),
        blastn: p,
        prodigal: None,
    }
}

fn run(
    delta: &str,
    recs: &[FastaRecord],
    genes: &[Gene],
    min_gained: u64,
) -> Vec<straincompass_types::GainedRow> {
    let d = parse_delta_str(delta).unwrap();
    let dir = std::env::temp_dir().join(format!("gained_test_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let (rows, _) = gained_regions(
        &no_tools(),
        &d,
        recs,
        genes,
        &dir.join("gained_regions.fa"),
        &dir,
        min_gained,
        false,
    )
    .unwrap();
    rows
}

#[test]
fn aligned_qry_intervals_merges_per_contig() {
    // Two overlapping blocks on ctg1 (the second reverse on the reference)
    // plus one on ctg2. The reverse block must still contribute ascending
    // query coordinates.
    let mut d = String::from("/ref.fa /qry.fa\nNUCMER\n");
    d.push_str(">chr1 ctg1 10000 5000\n1 1000 1 1000 0 0 0\n0\n");
    d.push_str(">chr1 ctg1 10000 5000\n3000 2001 900 1900 0 0 0\n0\n");
    d.push_str(">chr1 ctg2 10000 2000\n5000 5500 100 600 0 0 0\n0\n");
    let parsed = parse_delta_str(&d).unwrap();
    let ivs = parsed.aligned_qry_intervals();
    assert_eq!(ivs["ctg1"], vec![(1, 1900)], "overlapping blocks merge");
    assert_eq!(ivs["ctg2"], vec![(100, 600)]);
}

#[test]
fn gained_regions_complement_and_threshold() {
    let d = block(
        "chr1",
        "ctg1",
        10000,
        5000,
        &[(1, 1000, 1, 1000), (5001, 6000, 4001, 5000)],
    );
    let rows = run(&d, &[flat("ctg1", 5000)], &[], 500);
    assert_eq!(rows.len(), 1);
    assert_eq!(
        (rows[0].start, rows[0].end, rows[0].length),
        (1001, 4000, 3000)
    );
    assert!(!rows[0].at_contig_end, "flanked on both sides");

    let rows = run(&d, &[flat("ctg1", 5000)], &[], 5000);
    assert!(
        rows.is_empty(),
        "threshold above the region length drops it"
    );
}

#[test]
fn anchor_between_forward_flanks() {
    // Both flanks forward on chr1, ending at 1000 and starting at 1001:
    // a clean insertion point.
    let d = block(
        "chr1",
        "ctg1",
        10000,
        5000,
        &[(1, 1000, 1, 1000), (1001, 2000, 4001, 5000)],
    );
    let genes = vec![
        gene("chr1", "G1", 500, 1000),
        gene("chr1", "G2", 1001, 1500),
    ];
    let rows = run(&d, &[flat("ctg1", 5000)], &genes, 500);
    assert_eq!(rows.len(), 1);
    let r = &rows[0];
    assert_eq!(r.anchor, GainedAnchor::Between);
    assert_eq!(r.anchor_seqid, "chr1");
    assert_eq!((r.anchor_start, r.anchor_end), (1000, 1001));
    assert!(!r.flanks_disagree);
    assert_eq!(r.left_gene, "G1");
    assert_eq!(r.right_gene, "G2");
}

#[test]
fn anchor_reverse_flank_uses_the_low_reference_end() {
    // The left flank is reverse: query 1..1000 maps to reference 2000..1001
    // descending, so the query base just before the region (1000) sits at
    // reference 1001 - the block's LOW end, not its high one.
    let d = block("chr1", "ctg1", 10000, 5000, &[(1001, 2000, 1000, 1)]);
    let rows = run(&d, &[flat("ctg1", 5000)], &[], 500);
    assert_eq!(rows.len(), 1);
    let r = &rows[0];
    assert_eq!(r.anchor, GainedAnchor::Flank);
    assert_eq!(
        r.anchor_start, 1001,
        "a reverse left flank contributes ref_start, not ref_end"
    );
}

#[test]
fn anchor_flank_at_contig_end() {
    // The only block starts at query 2001, so query 1..2000 is gained and
    // has no left flank at all.
    let d = block("chr1", "ctg1", 10000, 5000, &[(5000, 7999, 2001, 5000)]);
    let rows = run(&d, &[flat("ctg1", 5000)], &[], 500);
    assert_eq!(rows.len(), 1);
    let r = &rows[0];
    assert_eq!((r.start, r.end), (1, 2000));
    assert_eq!(r.anchor, GainedAnchor::Flank);
    assert_eq!(
        r.anchor_start, 5000,
        "placed at the junction with its one flank"
    );
    assert!(r.at_contig_end);
}

#[test]
fn anchor_unanchored_whole_contig() {
    // ctg2 has no alignment anywhere: a plasmid, unplaceable.
    let d = block("chr1", "ctg1", 10000, 5000, &[(1, 5000, 1, 5000)]);
    let rows = run(&d, &[flat("ctg1", 5000), flat("ctg2", 3000)], &[], 500);
    assert_eq!(rows.len(), 1);
    let r = &rows[0];
    assert_eq!(r.qry_seqid, "ctg2");
    assert_eq!(r.anchor, GainedAnchor::Unanchored);
    assert_eq!(r.anchor_seqid, "");
    assert_eq!((r.anchor_start, r.anchor_end), (0, 0));
    assert!(r.at_contig_end, "the whole contig is the region");
}

#[test]
fn anchor_discordant_flanks_fall_back_to_the_left() {
    // The flanks land on different reference sequences, so no single
    // "between" position is honest.
    let mut d = String::from("/ref.fa /qry.fa\nNUCMER\n");
    d.push_str(">chr1 ctg1 10000 5000\n1 1000 1 1000 0 0 0\n0\n");
    d.push_str(">chr2 ctg1 10000 5000\n1 1000 4001 5000 0 0 0\n0\n");
    let rows = run(&d, &[flat("ctg1", 5000)], &[], 500);
    assert_eq!(rows.len(), 1);
    let r = &rows[0];
    assert_eq!(r.anchor, GainedAnchor::Flank);
    assert_eq!(r.anchor_seqid, "chr1", "falls back to the left flank");
    assert!(r.flanks_disagree);
}

#[test]
fn anchor_far_apart_flanks_disagree() {
    // Same reference sequence, but the flanks are 50 kb apart - further
    // than one insertion point can explain.
    let d = block(
        "chr1",
        "ctg1",
        200_000,
        5000,
        &[(1, 1000, 1, 1000), (51_001, 52_000, 4001, 5000)],
    );
    let rows = run(&d, &[flat("ctg1", 5000)], &[], 500);
    assert_eq!(rows[0].anchor, GainedAnchor::Flank);
    assert!(rows[0].flanks_disagree);
}

#[test]
fn gc_pct_ignores_ambiguous_bases() {
    // A 1000bp contig: the first 400 align, the rest is GC plus a run of N.
    let seq = format!("{}{}{}", "A".repeat(400), "GC".repeat(250), "N".repeat(100));
    let d = block("chr1", "ctg1", 10000, 1000, &[(1, 400, 1, 400)]);
    let rows = run(&d, &[rec("ctg1", &seq)], &[], 100);
    assert_eq!(rows.len(), 1);
    assert!(
        (rows[0].gc_pct - 100.0).abs() < 1e-9,
        "N is excluded from the denominator, not counted as not-GC: got {}",
        rows[0].gc_pct
    );
}

#[test]
fn gained_without_prodigal_degrades() {
    let d = block("chr1", "ctg1", 10000, 5000, &[(1, 1000, 1, 1000)]);
    let parsed = parse_delta_str(&d).unwrap();
    let dir = std::env::temp_dir().join("gained_test_degrade");
    std::fs::create_dir_all(&dir).unwrap();
    let (rows, status) = gained_regions(
        &no_tools(),
        &parsed,
        &[flat("ctg1", 5000)],
        &[],
        &dir.join("gained_regions.fa"),
        &dir,
        500,
        true,
    )
    .unwrap();
    assert!(!rows.is_empty());
    for r in &rows {
        assert_eq!(r.n_orfs, None, "no count is not a zero count");
        assert_eq!(r.n_orfs_complete, None);
        assert!(r.orfs.is_empty());
    }
    assert!(matches!(status, GainedOrfStatus::Unavailable(_)));
}

/// A stand-in gene finder: writes a fixed prodigal-style GFF to whatever
/// path follows `-o`. Real prodigal is not installed everywhere (and is
/// optional by design), but the invoke-and-parse path still has to be
/// exercised, including the coordinate shift and the partial flag.
fn stub_prodigal(dir: &std::path::Path, gff_body: &str) -> PathBuf {
    let gff = dir.join("stub_output.gff");
    std::fs::write(&gff, gff_body).unwrap();
    let script = dir.join("stub_prodigal.sh");
    std::fs::write(
        &script,
        format!(
            "#!/bin/sh\nwhile [ $# -gt 0 ]; do\n  if [ \"$1\" = \"-o\" ]; then cp '{}' \"$2\"; fi\n  shift\ndone\nexit 0\n",
            gff.display()
        ),
    )
    .unwrap();
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        std::fs::set_permissions(&script, std::fs::Permissions::from_mode(0o755)).unwrap();
    }
    script
}

#[test]
fn orfs_are_shifted_into_query_contig_coordinates() {
    let dir = std::env::temp_dir().join("gained_test_orfs");
    std::fs::create_dir_all(&dir).unwrap();
    // Region gr0 is query 1001..4000, so a gene at region offset 101..700
    // sits at query 1101..1700. The second call runs off an edge.
    let body = "##gff-version  3\n\
        # Sequence Data: seqnum=1;seqlen=3000\n\
        gr0\tProdigal_v2.6.3\tCDS\t101\t700\t45.6\t+\t0\tID=1_1;partial=00;conf=99.99\n\
        gr0\tProdigal_v2.6.3\tCDS\t2800\t3000\t12.1\t-\t0\tID=1_2;partial=01;conf=71.20\n";
    let prodigal = stub_prodigal(&dir, body);

    let d = block(
        "chr1",
        "ctg1",
        10000,
        5000,
        &[(1, 1000, 1, 1000), (5001, 6000, 4001, 5000)],
    );
    let parsed = parse_delta_str(&d).unwrap();
    let mut tools = no_tools();
    tools.prodigal = Some(prodigal);
    let (rows, status) = gained_regions(
        &tools,
        &parsed,
        &[flat("ctg1", 5000)],
        &[],
        &dir.join("gained_regions.fa"),
        &dir,
        500,
        true,
    )
    .unwrap();

    assert_eq!(status, GainedOrfStatus::Predicted);
    assert_eq!(rows.len(), 1);
    let r = &rows[0];
    assert_eq!((r.start, r.end), (1001, 4000));
    assert_eq!(r.n_orfs, Some(2));
    assert_eq!(
        r.n_orfs_complete,
        Some(1),
        "the partial call is not headlined"
    );
    assert_eq!((r.orfs[0].start, r.orfs[0].end), (1101, 1700));
    assert_eq!(r.orfs[0].strand, 1);
    assert!(!r.orfs[0].partial);
    assert!((r.orfs[0].confidence - 99.99).abs() < 1e-9);
    assert_eq!((r.orfs[1].start, r.orfs[1].end), (3800, 4000));
    assert_eq!(r.orfs[1].strand, -1);
    assert!(r.orfs[1].partial);

    // The region fasta is written for the gene finder and kept for the user.
    let fa = std::fs::read_to_string(dir.join("gained_regions.fa")).unwrap();
    assert!(
        fa.starts_with(">gr0 ctg1:1001-4000\n"),
        "got {:?}",
        &fa[..40.min(fa.len())]
    );
}

#[test]
fn a_failing_gene_finder_does_not_fail_the_comparison() {
    let dir = std::env::temp_dir().join("gained_test_fail");
    std::fs::create_dir_all(&dir).unwrap();
    let script = dir.join("broken.sh");
    std::fs::write(&script, "#!/bin/sh\necho 'boom' >&2\nexit 1\n").unwrap();
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        std::fs::set_permissions(&script, std::fs::Permissions::from_mode(0o755)).unwrap();
    }
    let d = block("chr1", "ctg1", 10000, 5000, &[(1, 1000, 1, 1000)]);
    let parsed = parse_delta_str(&d).unwrap();
    let mut tools = no_tools();
    tools.prodigal = Some(script);
    let (rows, status) = gained_regions(
        &tools,
        &parsed,
        &[flat("ctg1", 5000)],
        &[],
        &dir.join("gained_regions.fa"),
        &dir,
        500,
        true,
    )
    .unwrap();
    assert!(
        !rows.is_empty(),
        "the regions survive a gene finder failure"
    );
    assert_eq!(rows[0].n_orfs, None);
    match status {
        GainedOrfStatus::Unavailable(m) => assert!(m.contains("boom"), "got {m}"),
        other => panic!("expected Unavailable, got {other:?}"),
    }
}
