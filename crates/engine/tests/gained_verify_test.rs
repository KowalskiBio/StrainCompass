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

/// A stand-in blastn: writes fixed outfmt-6 lines to whatever path
/// follows -out.
fn stub_blastn(dir: &Path, hits: &str) -> PathBuf {
    let p = dir.join("stub_blastn.sh");
    write_executable(
        &p,
        &format!(
            "#!/bin/sh\nout=\"\"\nprev=\"\"\nfor a in \"$@\"; do\n  if [ \"$prev\" = \"-out\" ]; then out=\"$a\"; fi\n  prev=\"$a\"\ndone\ncat > \"$out\" <<'HITS'\n{hits}HITS\n"
        ),
    );
    p
}

fn check(dir: &Path, hits: &str) -> straincompass_types::GainedVerify {
    std::fs::write(dir.join("ref.fa"), ">chr1\nACGTACGT\n>chr2\nTTTTTTTT\n").unwrap();
    std::fs::write(dir.join("qry.fa"), ">ctg1\nACGTACGTAC\n").unwrap();
    let mut t = tools(dir);
    t.blastn = stub_blastn(dir, hits);
    gained_verify(
        &t,
        &dir.join("ref.fa"),
        &dir.join("qry.fa"),
        "ctg1",
        3,
        8,
        &dir.join("work"),
    )
    .unwrap()
}

fn tools(dir: &Path) -> ToolPaths {
    let n = PathBuf::from("/nonexistent");
    ToolPaths {
        nucmer: n.clone(),
        show_coords: n.clone(),
        show_snps: n.clone(),
        dnadiff: n.clone(),
        makeblastdb: stub_makeblastdb(dir),
        blastn: stub_blastn(dir, ""),
        prodigal: None,
    }
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
fn no_hits_is_an_empty_table_not_an_error() {
    let d = test_dir("empty");
    let v = check(&d, "");
    assert!(v.hits.is_empty());
    let _ = std::fs::remove_dir_all(&d);
}

#[test]
fn an_unknown_contig_is_a_friendly_error() {
    let d = test_dir("unknown");
    std::fs::write(d.join("ref.fa"), ">chr1\nACGTACGT\n").unwrap();
    std::fs::write(d.join("qry.fa"), ">ctg1\nACGTACGTAC\n").unwrap();
    let err = gained_verify(
        &tools(&d),
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
