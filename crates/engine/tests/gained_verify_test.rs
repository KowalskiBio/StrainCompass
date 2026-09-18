//! The reference back-check of gained regions, exercised through stand-in
//! blast binaries: real BLAST+ is not installed everywhere (and is not on
//! the CI runner), the same trade the prodigal stubs in gained_test make.

use std::path::{Path, PathBuf};
use straincompass_engine::blast::gained_verify;
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
/// ceiling the caller passed (field 8 of outfmt 6) and writes the
/// survivors to whatever path follows -out. Filtering matters: the
/// engine picks the ceilings, and without it the stub would report hits
/// the real tool would have withheld.
fn stub_hits_script(dir: &Path, name: &str, hits: &str) -> PathBuf {
    let p = dir.join(format!("stub_{name}.sh"));
    write_executable(
        &p,
        &format!(
            "#!/bin/sh\nout=\"\"\nev=\"\"\nprev=\"\"\nfor a in \"$@\"; do\n  if [ \"$prev\" = \"-out\" ]; then out=\"$a\"; fi\n  if [ \"$prev\" = \"-evalue\" ]; then ev=\"$a\"; fi\n  prev=\"$a\"\ndone\nawk -v ev=\"$ev\" 'NF == 0 || $8+0 <= ev+0' > \"$out\" <<'HITS'\n{hits}HITS\n"
        ),
    );
    p
}

/// One back-check with the given blastn and tblastx hit tables. The
/// reference is chr1/chr2, the query contig ctg1 ("ACGTACGTAC"), and the
/// region under test is always ctg1:3-8 ("GTACGT").
fn check(dir: &Path, nuc_hits: &str, tx_hits: &str) -> straincompass_types::GainedVerify {
    std::fs::write(dir.join("ref.fa"), ">chr1\nACGTACGT\n>chr2\nTTTTTTTT\n").unwrap();
    std::fs::write(dir.join("qry.fa"), ">ctg1\nACGTACGTAC\n").unwrap();
    let tools = ToolPaths {
        nucmer: PathBuf::from("/nonexistent"),
        show_coords: PathBuf::from("/nonexistent"),
        show_snps: PathBuf::from("/nonexistent"),
        dnadiff: PathBuf::from("/nonexistent"),
        makeblastdb: stub_makeblastdb(dir),
        blastn: stub_hits_script(dir, "blastn", nuc_hits),
        tblastx: stub_hits_script(dir, "tblastx", tx_hits),
        prodigal: None,
    };
    gained_verify(
        &tools,
        &dir.join("ref.fa"),
        &dir.join("qry.fa"),
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
        blastn: stub_hits_script(&d, "blastn", ""),
        tblastx: stub_hits_script(&d, "tblastx", "chr1\t5\t50\t45.0\t15\t1\t45\t1e-9\t60\n"),
        prodigal: None,
    };
    let v = gained_verify(
        &tools,
        &d.join("ref.fa"),
        &d.join("qry.fa"),
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
        blastn: stub_hits_script(&d, "blastn", ""),
        tblastx: stub_hits_script(&d, "tblastx", ""),
        prodigal: None,
    };
    let err = gained_verify(
        &tools,
        &d.join("ref.fa"),
        &d.join("qry.fa"),
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
