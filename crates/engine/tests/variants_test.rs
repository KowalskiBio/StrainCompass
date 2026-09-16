//! Variant event extraction tests (pure, no external tools).
//!
//! Delta semantics used by the fixtures (see delta.rs): each delta value
//! is one gap column; `|d| - 1` paired bases come first. `d < 0` = one
//! inserted query base (gap in the reference row), `d > 0` = one deleted
//! reference base (gap in the query row).

use bactiment_engine::delta::parse_delta_str;
use bactiment_engine::fasta::parse_fasta_str;
use bactiment_engine::variants::variant_events;

const REF_FASTA: &str = ">chr1\nACGTACGT\n";

fn delta(body: &str) -> String {
    format!(">chr1 q1\n{body}\n")
}

fn events(text: &str) -> std::collections::BTreeMap<String, bactiment_types::AlignmentEvents> {
    let delta = parse_delta_str(text).unwrap();
    let ref_records = parse_fasta_str(REF_FASTA).unwrap();
    let qry_records = parse_fasta_str(">q1\nGGAAGT\n").unwrap();
    variant_events(&delta, &ref_records, &qry_records).unwrap()
}

#[test]
fn forward_block_snps_dels_and_ins() {
    // Block: ref 1..8, qry 1..6, deltas [-1, 1, 1, 2]:
    //   -1      -> insertion of qry base 1 before ref base 1 (pos 0)
    //   1, 1    -> deletion of ref bases 1-2 (one merged event)
    //   2       -> 1 paired (ref3/qry2), then deletion of ref 4
    //   rest    -> paired ref5..8 with qry3..6
    // qry "GGAAGT": G inserted, G matches ref3, A matches ref5,
    // A mismatches ref6 (SNP C->A), G matches ref7, T matches ref8.
    let d = delta("1 8 1 6 1 0 0\n-1\n1\n1\n2\n0");
    let ev = events(&d);
    let e = &ev["chr1"];
    assert_eq!(e.ins.len(), 1);
    assert_eq!(e.ins[0].pos, 0);
    assert_eq!(e.ins[0].seq, "G");
    assert_eq!(e.dels.len(), 2);
    assert_eq!(e.dels[0].pos, 1);
    assert_eq!(e.dels[0].len, 2);
    assert_eq!(e.dels[1].pos, 4);
    assert_eq!(e.dels[1].len, 1);
    assert_eq!(e.snps.len(), 1);
    assert_eq!(e.snps[0].pos, 6);
    assert_eq!(e.snps[0].r, b'C');
    assert_eq!(e.snps[0].q, b'A');
}

#[test]
fn consecutive_insertions_merge_into_one_event() {
    // deltas [-1, -1, -1]: three inserted query bases before ref base 1,
    // then the remaining 8 paired bases all match.
    let d = delta("1 8 1 11 0 0 0\n-1\n-1\n-1\n0");
    let delta = parse_delta_str(&d).unwrap();
    let ref_records = parse_fasta_str(REF_FASTA).unwrap();
    // qry = 3 inserted + 8 paired = 11 bases; the paired bases match.
    let qry_records = parse_fasta_str(">q1\nTTAACGTACGT\n").unwrap();
    let ev = variant_events(&delta, &ref_records, &qry_records).unwrap();
    let e = &ev["chr1"];
    assert_eq!(e.ins.len(), 1);
    assert_eq!(e.ins[0].pos, 0);
    assert_eq!(e.ins[0].seq, "TTA");
    assert!(e.snps.is_empty(), "{:?}", e.snps);
    assert!(e.dels.is_empty(), "{:?}", e.dels);
}

#[test]
fn reverse_block_reports_query_bases_in_reference_orientation() {
    // Same geometry as the forward test but on the reverse strand
    // (header s2 > e2, qry 1..6). The forward query sequence
    // "ACTTCC" reverse complements to "GGAAGT", so the block must
    // produce exactly the same display events as the forward test
    // with a forward query "GGAAGT".
    let text = ">chr1 q2\n1 8 6 1 1 0 0\n-1\n1\n1\n2\n0\n";
    let delta = parse_delta_str(text).unwrap();
    let ref_records = parse_fasta_str(REF_FASTA).unwrap();
    let qry_records = parse_fasta_str(">q2\nACTTCC\n").unwrap();
    let ev = variant_events(&delta, &ref_records, &qry_records).unwrap();
    let e = &ev["chr1"];
    assert_eq!(e.ins.len(), 1);
    assert_eq!(e.ins[0].pos, 0);
    // the insertion sits at the block start in display space, i.e. at
    // the block end of the forward query: forward 'C' displays as 'G'
    assert_eq!(e.ins[0].seq, "G");
    assert_eq!(e.dels.len(), 2);
    assert_eq!((e.dels[0].pos, e.dels[0].len), (1, 2));
    assert_eq!((e.dels[1].pos, e.dels[1].len), (4, 1));
    assert_eq!(e.snps.len(), 1);
    assert_eq!(e.snps[0].pos, 6);
    assert_eq!((e.snps[0].r, e.snps[0].q), (b'C', b'A'));
}

#[test]
fn overlapping_blocks_report_each_event_once() {
    // Two identical blocks covering the same reference stretch: the
    // events must not be doubled.
    let text = ">chr1 q1\n1 8 1 6 1 0 0\n-1\n1\n1\n2\n0\n1 8 1 6 1 0 0\n-1\n1\n1\n2\n0\n";
    let ev = events(text);
    let e = &ev["chr1"];
    assert_eq!(e.ins.len(), 1);
    assert_eq!(e.dels.len(), 2);
    assert_eq!(e.snps.len(), 1);
}

#[test]
fn blocks_on_different_reference_seqids_are_kept_apart() {
    // chrA: 4 paired mismatches (AAAA vs CCCC). chrB: one deletion
    // from delta 3 (2 paired, then ref base 3 deleted).
    let text = ">chrA q1\n1 4 1 4 4 0 0\n0\n>chrB q1\n1 4 1 4 1 0 0\n3\n0\n";
    let delta = parse_delta_str(text).unwrap();
    let ref_records = parse_fasta_str(">chrA\nAAAA\n>chrB\nCCCC\n").unwrap();
    let qry_records = parse_fasta_str(">q1\nCCCC\n").unwrap();
    let ev = variant_events(&delta, &ref_records, &qry_records).unwrap();
    assert!(ev.contains_key("chrA"));
    assert!(ev.contains_key("chrB"));
    assert_eq!(ev["chrA"].snps.len(), 4);
    assert_eq!(ev["chrA"].dels.len(), 0);
    assert_eq!(ev["chrB"].snps.len(), 0);
    assert_eq!(ev["chrB"].dels.len(), 1);
    assert_eq!((ev["chrB"].dels[0].pos, ev["chrB"].dels[0].len), (3, 1));
}
