use bactiment_engine::delta::parse_delta_str;

/// The exact delta output of a controlled nucmer run:
/// reference 3000bp, query 3003bp with 5 SNPs (ref 301,601,901,1201,1501),
/// a 10bp insertion in the query after ref position 1800, and a 7bp
/// deletion in the query covering ref 2391-2397.
/// show-coords reported 99.27% identity for this alignment.
const DELTA: &str = "/ref.fa /qry.fa
NUCMER
>chr1 ctg1 3000 3003
1 3000 1 3003 22 22 0
-1801
-1
-1
-1
-1
-1
-1
-1
-1
-1
591
1
1
1
1
1
1
0
";

#[test]
fn parses_and_scores_like_show_coords() {
    let d = parse_delta_str(DELTA).unwrap();
    assert_eq!(d.alignments.len(), 1);
    let a = &d.alignments[0];
    assert_eq!(a.ref_seqid, "chr1");
    assert_eq!(a.qry_seqid, "ctg1");
    assert_eq!(a.ref_start, 1);
    assert_eq!(a.ref_end, 3000);
    assert_eq!(a.qry_lo, 1);
    assert_eq!(a.qry_hi, 3003);
    assert!(!a.qry_rev);
    assert_eq!(a.errors, 22);
    assert_eq!(a.mismatches(), 5);
    assert_eq!(a.qry_ins(), 10);
    assert_eq!(a.ref_ins(), 7);
    assert_eq!(a.paired(), 2993);
    assert_eq!(a.aligned_len(), 3010);
    let idy = a.identity();
    assert!((idy - 99.269).abs() < 0.01, "identity {idy}");
}

#[test]
fn reconstructs_alignment_rows() {
    use bactiment_engine::msa;
    let d = parse_delta_str(DELTA).unwrap();
    let a = &d.alignments[0];
    // Deterministic reference and the matching query: build ref, then apply
    // the same edits used to produce this delta.
    let ref_seq: Vec<u8> = (0..3000)
        .map(|i| b"ACGT"[i % 4])
        .collect();
    // query = ref with SNPs at 0-based 300,600,900,1200,1500, insertion of
    // 10 bases after index 1800, deletion of ref bases 2390..2397 (0-based).
    let mut q: Vec<u8> = ref_seq.clone();
    for &p in &[300usize, 600, 900, 1200, 1500] {
        let orig = q[p];
        q[p] = match orig {
            b'A' => b'C',
            b'C' => b'G',
            b'G' => b'T',
            _ => b'A',
        };
    }
    let insert = vec![b'A'; 10];
    let mut q: Vec<u8> = q[..1800].to_vec();
    q.extend_from_slice(&insert);
    q.extend_from_slice(&ref_seq[1800..2390]);
    q.extend_from_slice(&ref_seq[2397..]);
    assert_eq!(q.len(), 3003);

    let pw = msa::reconstruct(a, &ref_seq, &q);
    assert_eq!(pw.len(), 3010);
    // ref row: 3000 bases + 10 dashes
    let dashes_ref = pw.ref_row.iter().filter(|&&c| c == b'-').count();
    let dashes_qry = pw.qry_row.iter().filter(|&&c| c == b'-').count();
    assert_eq!(dashes_ref, 10);
    assert_eq!(dashes_qry, 7);
    // mismatch count
    let mm = pw
        .ref_row
        .iter()
        .zip(&pw.qry_row)
        .filter(|(r, s)| **r != b'-' && **s != b'-' && r != s)
        .count();
    assert_eq!(mm, 5);
    // the gap in the reference row sits right after ref position 1800
    let idx = pw
        .ref_row
        .windows(10)
        .position(|w| w.iter().all(|&c| c == b'-'))
        .expect("10 dash run in ref row");
    assert_eq!(pw.ref_pos[idx], 1800);
    // query positions of the inserted run are 1801..1810
    assert_eq!(pw.qry_pos[idx], 1801);
    assert_eq!(pw.qry_pos[idx + 9], 1810);
    // the gap in the query row starts after ref 2390 (deleted ref 2391-2397)
    let qidx = pw
        .qry_row
        .windows(7)
        .position(|w| w.iter().all(|&c| c == b'-'))
        .expect("7 dash run in qry row");
    assert_eq!(pw.qry_row[qidx..qidx + 7].iter().all(|&c| c == b'-'), true);
    assert_eq!(pw.ref_pos[qidx], 2391);
    assert_eq!(pw.ref_pos[qidx + 6], 2397);
}

#[test]
fn merged_intervals() {
    let mut ivs = vec![(10u64, 20u64), (5, 8), (21, 30), (40, 50)];
    bactiment_engine::delta::merge_intervals(&mut ivs);
    assert_eq!(ivs, vec![(5, 8), (10, 30), (40, 50)]);
}

#[test]
fn reverse_strand_header() {
    // s2 > e2 means reverse alignment
    let d = parse_delta_str(">r q 100 100\n1 100 100 1 0 0 0\n0\n").unwrap();
    assert!(d.alignments[0].qry_rev);
    assert_eq!(d.alignments[0].qry_lo, 1);
    assert_eq!(d.alignments[0].qry_hi, 100);
}
