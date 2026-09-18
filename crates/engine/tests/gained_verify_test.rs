//! The reference back-check of gained regions, exercised through stand-in
//! blast binaries: real BLAST+ is not installed everywhere (and is not on
//! the CI runner), the same trade the prodigal stubs in gained_test make.

use std::path::{Path, PathBuf};
use straincompass_engine::blast::{gained_verify, RegionInputs};
use straincompass_engine::tools::ToolPaths;

fn write_executable(path: &Path, body: &str) {
    std::fs::write(path, body).unwrap();
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        std::fs::set_permissions(path, std::fs::Permissions::from_mode(0o755)).unwrap();
    }
}

/// A stand-in makeblastdb: the back-check never reads the database it
/// builds, only whether the build succeeded.
fn stub_makeblastdb(dir: &Path) -> PathBuf {
    let p = dir.join("stub_makeblastdb.sh");
    write_executable(&p, "#!/bin/sh\nexit 0\n");
    p
}

/// A stand-in blast tool: filters its fixed hit table to the E-value
/// ceiling the caller passed (`evfield` is the column of outfmt 6 the
/// engine asked the real tool for, different per outfmt) and writes the
/// survivors to whatever path follows -out. Filtering matters: the
/// engine picks the ceilings, and without it the stub would report hits
/// the real tool would have withheld.
fn stub_hits_script(dir: &Path, name: &str, hits: &str, evfield: usize) -> PathBuf {
    let p = dir.join(format!("stub_{name}.sh"));
    write_executable(
        &p,
        &format!(
            "#!/bin/sh\nout=\"\"\nev=\"\"\nprev=\"\"\nfor a in \"$@\"; do\n  if [ \"$prev\" = \"-out\" ]; then out=\"$a\"; fi\n  if [ \"$prev\" = \"-evalue\" ]; then ev=\"$a\"; fi\n  prev=\"$a\"\ndone\nawk -v ev=\"$ev\" -v ef={evfield} 'NF == 0 || $ef+0 <= ev+0' > \"$out\" <<'HITS'\n{hits}HITS\n"
        ),
    );
    p
}

/// One back-check with the given blastn and tblastx hit tables. The
/// reference is chr1/chr2, the query contig ctg1 ("ACGTACGTAC"), and the
/// region under test is always ctg1:3-8 ("GTACGT"). The reference is
/// annotated with two genes so the hit labeling is exercised: G1 with a
/// symbol covers chr1:100-200, G2 without one covers chr2:700-950.
fn check(dir: &Path, nuc_hits: &str, tx_hits: &str) -> straincompass_types::GainedVerify {
    std::fs::write(dir.join("ref.fa"), ">chr1\nACGTACGT\n>chr2\nTTTTTTTT\n").unwrap();
    std::fs::write(
        dir.join("ref.gff"),
        "##gff-version 3\n\
         chr1\t.\tgene\t100\t200\t.\t+\t.\tlocus_tag=G1;gene=glx;gene_biotype=protein_coding\n\
         chr1\t.\tCDS\t100\t200\t.\t+\t0\tlocus_tag=G1;product=glucose oxidase;protein_id=WP_G1\n\
         chr2\t.\tgene\t700\t950\t.\t-\t.\tlocus_tag=G2;gene_biotype=protein_coding\n\
         chr2\t.\tCDS\t700\t950\t.\t-\t0\tlocus_tag=G2;product=minor pseudopilin;protein_id=WP_G2\n",
    )
    .unwrap();
    std::fs::write(dir.join("qry.fa"), ">ctg1\nACGTACGTAC\n").unwrap();
    let tools = ToolPaths {
        nucmer: PathBuf::from("/nonexistent"),
        show_coords: PathBuf::from("/nonexistent"),
        show_snps: PathBuf::from("/nonexistent"),
        dnadiff: PathBuf::from("/nonexistent"),
        makeblastdb: stub_makeblastdb(dir),
        blastn: stub_hits_script(dir, "blastn", nuc_hits, 8),
        tblastx: stub_hits_script(dir, "tblastx", tx_hits, 8),
        blastx: stub_hits_script(dir, "blastx", "", 5),
        prodigal: None,
    };
    gained_verify(
        &tools,
        &RegionInputs {
            ref_fasta: &dir.join("ref.fa"),
            ref_gff: &dir.join("ref.gff"),
            qry_fasta: &dir.join("qry.fa"),
        },
        "ctg1",
        3,
        8,
        &dir.join("work"),
    )
    .unwrap()
}

fn test_dir(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("gained_verify_{name}_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

#[test]
fn hits_are_parsed_sorted_and_normalized() {
    let d = test_dir("parsed");
    // The weaker hit is listed first and its subject span runs backwards
    // (sstart > send): it must come second and carry an ordered interval.
    let v = check(
        &d,
        "chr2\t900\t750\t98.1\t151\t50\t200\t2e-20\t120\nchr1\t100\t200\t99.5\t101\t1\t101\t1e-30\t185\n",
        "",
    );
    assert_eq!((v.start, v.end), (3, 8));
    assert_eq!(v.hits.len(), 2);
    assert_eq!(v.hits[0].ref_seqid, "chr1", "best bitscore first");
    assert_eq!((v.hits[0].ref_start, v.hits[0].ref_end), (100, 200));
    assert_eq!(v.hits[0].identity, 99.5);
    assert_eq!(
        (v.hits[1].ref_start, v.hits[1].ref_end),
        (750, 900),
        "a reverse-strand hit arrives as sstart > send and must be ordered"
    );
    // Every hit says which genes it overlaps: the symbol when the
    // annotation has one, the product when it does not.
    assert_eq!(v.hits[0].genes, vec!["G1 (glx)"]);
    assert_eq!(v.hits[1].genes, vec!["G2 (minor pseudopilin)"]);
    let _ = std::fs::remove_dir_all(&d);
}

#[test]
fn the_evalue_ceiling_splits_the_tiers() {
    let d = test_dir("tiers");
    // 2e-20 and 1e-30 are at or below 1e-5 (strong); 0.5 and 3 are
    // above it but inside the loose ceiling of 10 (weak); 500 is outside
    // the search entirely and must not be reported at all - the stub
    // would only emit it if the ceiling were wrong, so asserting its
    // absence pins the parameter.
    let v = check(
        &d,
        "chr1\t10\t20\t99.0\t11\t1\t11\t2e-20\t120\n\
         chr1\t30\t40\t98.0\t11\t1\t11\t1e-30\t185\n\
         chr1\t50\t60\t90.0\t11\t1\t11\t0.5\t25\n\
         chr1\t70\t80\t85.0\t11\t1\t11\t3\t20\n\
         chr1\t90\t100\t80.0\t11\t1\t11\t500\t10\n",
        "",
    );
    assert_eq!(v.hits.len(), 2, "strong tier keeps E <= 1e-5");
    assert!(v.hits.iter().all(|h| h.evalue <= 1e-5));
    assert_eq!(v.weak_hits.len(), 2, "weak tier keeps 1e-5 < E <= 10");
    assert!(v.weak_hits.iter().all(|h| h.evalue > 1e-5));
    let _ = std::fs::remove_dir_all(&d);
}

#[test]
fn no_hits_is_an_empty_table_not_an_error() {
    let d = test_dir("empty");
    let v = check(&d, "", "");
    assert!(v.hits.is_empty());
    assert!(v.weak_hits.is_empty());
    let _ = std::fs::remove_dir_all(&d);
}

#[test]
fn the_translated_search_seconds_a_silent_nucleotide_search() {
    let d = test_dir("tx");
    let v = check(&d, "", "chr2\t5\t50\t45.0\t15\t1\t45\t1e-9\t60\n");
    assert!(
        v.hits.is_empty(),
        "the nucleotide search is silent so the translated one runs"
    );
    let tx = v.tx_hits.expect("the translated search ran");
    assert_eq!(tx.len(), 1);
    assert_eq!(tx[0].ref_seqid, "chr2");
    assert_eq!(tx[0].identity, 45.0, "tblastx identity is over amino acids");
    let _ = std::fs::remove_dir_all(&d);
}

#[test]
fn the_translated_search_skips_a_region_with_strong_hits() {
    let d = test_dir("tx_skip");
    let v = check(
        &d,
        "chr1\t100\t200\t99.5\t101\t1\t101\t1e-30\t185\n",
        "chr2\t5\t50\t45.0\t15\t1\t45\t1e-9\t60\n",
    );
    assert_eq!(v.hits.len(), 1);
    assert!(
        v.tx_hits.is_none(),
        "nothing to second-guess: the strong tier already answered"
    );
    assert!(v.tx_note.is_none());
    let _ = std::fs::remove_dir_all(&d);
}

#[test]
fn an_overlong_region_skips_the_translated_search_with_a_note() {
    let d = test_dir("tx_long");
    // A region one base over the cap, so the skip is exercised exactly
    // at the boundary rather than far beyond it.
    let long = "A".repeat(50_001);
    std::fs::write(d.join("ref.fa"), ">chr1\nACGTACGT\n").unwrap();
    std::fs::write(d.join("qry.fa"), format!(">ctg1\n{long}\n")).unwrap();
    let tools = ToolPaths {
        nucmer: PathBuf::from("/nonexistent"),
        show_coords: PathBuf::from("/nonexistent"),
        show_snps: PathBuf::from("/nonexistent"),
        dnadiff: PathBuf::from("/nonexistent"),
        makeblastdb: stub_makeblastdb(&d),
        blastn: stub_hits_script(&d, "blastn", "", 8),
        tblastx: stub_hits_script(&d, "tblastx", "chr1\t5\t50\t45.0\t15\t1\t45\t1e-9\t60\n", 8),
        blastx: stub_hits_script(&d, "blastx", "", 5),
        prodigal: None,
    };
    let v = gained_verify(
        &tools,
        &RegionInputs {
            ref_fasta: &d.join("ref.fa"),
            ref_gff: &d.join("ref.gff"),
            qry_fasta: &d.join("qry.fa"),
        },
        "ctg1",
        1,
        50_001,
        &d.join("work"),
    )
    .unwrap();
    assert!(v.hits.is_empty());
    assert!(v.tx_hits.is_none());
    let note = v.tx_note.expect("the skip is explained, not silent");
    assert!(note.contains("50 kb"), "got {note}");
    let _ = std::fs::remove_dir_all(&d);
}

#[test]
fn the_longest_exact_match_is_reported_with_its_places() {
    let d = test_dir("lcs");
    // Region 3-8 of ctg1 is "GTACGT"; chr1 is "ACGTACGT" and contains
    // the whole of it, so the longest exact match is the region itself.
    let v = check(&d, "", "");
    assert_eq!(v.longest_exact_bp, 6);
    assert_eq!(v.longest_exact_seqid, "chr1");
    assert_eq!((v.longest_exact_start, v.longest_exact_end), (3, 8));
    assert_eq!(
        (v.longest_exact_qry_start, v.longest_exact_qry_end),
        (1, 6),
        "region-relative coordinates on the region side"
    );
    let _ = std::fs::remove_dir_all(&d);
}

#[test]
fn an_unknown_contig_is_a_friendly_error() {
    let d = test_dir("unknown");
    std::fs::write(d.join("ref.fa"), ">chr1\nACGTACGT\n").unwrap();
    std::fs::write(d.join("qry.fa"), ">ctg1\nACGTACGTAC\n").unwrap();
    let tools = ToolPaths {
        nucmer: PathBuf::from("/nonexistent"),
        show_coords: PathBuf::from("/nonexistent"),
        show_snps: PathBuf::from("/nonexistent"),
        dnadiff: PathBuf::from("/nonexistent"),
        makeblastdb: stub_makeblastdb(&d),
        blastn: stub_hits_script(&d, "blastn", "", 8),
        tblastx: stub_hits_script(&d, "tblastx", "", 8),
        blastx: stub_hits_script(&d, "blastx", "", 5),
        prodigal: None,
    };
    let err = gained_verify(
        &tools,
        &RegionInputs {
            ref_fasta: &d.join("ref.fa"),
            ref_gff: &d.join("ref.gff"),
            qry_fasta: &d.join("qry.fa"),
        },
        "nope",
        1,
        5,
        &d.join("work"),
    )
    .unwrap_err();
    let msg = format!("{err:?}");
    assert!(msg.contains("not in this query's fasta"), "got {msg}");
    let _ = std::fs::remove_dir_all(&d);
}

/// One identification run over the check() fixtures: reference chr1
/// ("ACGTACGT", gene G1 = protein "TY" from its two clean codons) and
/// query ctg1 ("ACGTACGTAC"), with one predicted ORF spanning ctg1:1-9.
fn identify(
    dir: &Path,
    orfs: &[straincompass_types::GainedOrf],
    blastx_hits: &str,
) -> straincompass_types::GainedIdentify {
    std::fs::write(dir.join("ref.fa"), ">chr1\nACGTACGT\n").unwrap();
    std::fs::write(
        dir.join("ref.gff"),
        "##gff-version 3\n\
         chr1\t.\tgene\t1\t8\t.\t+\t.\tlocus_tag=G1;gene=glx;gene_biotype=protein_coding\n\
         chr1\t.\tCDS\t1\t8\t.\t+\t0\tlocus_tag=G1;product=glucose oxidase;protein_id=WP_G1\n",
    )
    .unwrap();
    std::fs::write(dir.join("qry.fa"), ">ctg1\nACGTACGTAC\n").unwrap();
    let tools = ToolPaths {
        nucmer: PathBuf::from("/nonexistent"),
        show_coords: PathBuf::from("/nonexistent"),
        show_snps: PathBuf::from("/nonexistent"),
        dnadiff: PathBuf::from("/nonexistent"),
        makeblastdb: stub_makeblastdb(dir),
        blastn: stub_hits_script(dir, "blastn", "", 8),
        tblastx: stub_hits_script(dir, "tblastx", "", 8),
        blastx: stub_hits_script(dir, "blastx", blastx_hits, 5),
        prodigal: None,
    };
    straincompass_engine::blast::gained_identify(
        &tools,
        &RegionInputs {
            ref_fasta: &dir.join("ref.fa"),
            ref_gff: &dir.join("ref.gff"),
            qry_fasta: &dir.join("qry.fa"),
        },
        "ctg1",
        1,
        9,
        orfs,
        &dir.join("work"),
    )
    .unwrap()
}

#[test]
fn the_best_reference_protein_names_the_orf() {
    let d = test_dir("identify_hit");
    let orfs = vec![straincompass_types::GainedOrf {
        start: 1,
        end: 9,
        strand: 1,
        partial: false,
        confidence: 99.0,
    }];
    let v = identify(
        &d,
        &orfs,
        // Two hits for the same ORF: the better bitscore wins.
        "o0\tG1\t80.0\t70.0\t1e-5\t40\no0\tG1\t99.0\t100.0\t1e-30\t50\n",
    );
    assert_eq!(v.orfs.len(), 1);
    assert_eq!(v.orfs[0].seq, "ACGTACGTA");
    let m = v.orfs[0].best.as_ref().expect("the ORF must be named");
    assert_eq!(m.locus_tag, "G1");
    assert_eq!(m.protein_id, "WP_G1");
    assert_eq!(m.label, "glx");
    assert_eq!(m.identity, 99.0);
    assert_eq!(m.coverage, 100.0);
    assert_eq!(m.evalue, 1e-30);
    assert_eq!(v.region_seq, "ACGTACGTA");
    let _ = std::fs::remove_dir_all(&d);
}

#[test]
fn an_orf_without_a_reference_match_stays_unnamed() {
    let d = test_dir("identify_miss");
    let orfs = vec![
        straincompass_types::GainedOrf {
            start: 1,
            end: 9,
            strand: 1,
            partial: false,
            confidence: 90.0,
        },
        straincompass_types::GainedOrf {
            start: 2,
            end: 8,
            strand: -1,
            partial: true,
            confidence: 50.0,
        },
    ];
    let v = identify(&d, &orfs, "");
    assert_eq!(v.orfs.len(), 2);
    assert!(v.orfs[0].best.is_none());
    // A reverse-strand ORF still arrives as its own plus sequence.
    assert_eq!(v.orfs[1].seq, "ACGTACG");
    let _ = std::fs::remove_dir_all(&d);
}

#[test]
fn a_region_without_orfs_is_answered_without_searching() {
    let d = test_dir("identify_noorfs");
    // The blastx stub would fail the run if the engine called it, which
    // is the point: there is nothing to name, so nothing may run.
    std::fs::write(d.join("sentinel"), "").unwrap();
    let v = identify(&d, &[], "");
    assert!(v.orfs.is_empty());
    assert_eq!(v.region_seq, "ACGTACGTA");
    let _ = std::fs::remove_dir_all(&d);
}
