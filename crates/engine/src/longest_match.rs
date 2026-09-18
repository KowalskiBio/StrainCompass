//! The longest stretch of bases two sequences share exactly, found
//! without BLAST and without an E-value.
//!
//! The reference back-check of a gained region needs this because every
//! search has a sensitivity floor, while "what is the longest exact
//! match anywhere" is a fact about the two strings - reproducible,
//! cutoff-free, and interpretable against the ~log4(n*m) that chance
//! alone produces for unrelated DNA of those sizes.
//!
//! A suffix automaton is built over the *region* (the small side, at
//! most a few hundred kb), then every reference contig is scanned
//! through it once: O(region + reference) time and O(region) memory. An
//! automaton over the reference instead would be several times the
//! reference in states and, on a genome-sized input, needlessly heavy
//! for a request handler.

use crate::fasta::FastaRecord;

/// Where the longest exact match was found: on the reference and on the
/// region itself, both 1-based inclusive. Several matches may share the
/// maximum length; one example is reported.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LongestMatch {
    pub length: u64,
    pub ref_seqid: String,
    pub ref_start: u64,
    pub ref_end: u64,
    pub qry_start: u64,
    pub qry_end: u64,
}

/// The longest common substring of `needle` (the region) and any
/// position of the haystack records. When the two share nothing at all,
/// `length` is 0 and the positions are all 0.
pub fn longest_exact_match(needle: &[u8], haystack: &[FastaRecord]) -> LongestMatch {
    // Case-insensitive on purpose: FASTA case carries no meaning here,
    // and a soft-masked reference would otherwise hide its own matches.
    let needle: Vec<u8> = needle.iter().map(|b| b.to_ascii_uppercase()).collect();
    let sam = SuffixAutomaton::new(&needle);

    let mut best = LongestMatch {
        length: 0,
        ref_seqid: String::new(),
        ref_start: 0,
        ref_end: 0,
        qry_start: 0,
        qry_end: 0,
    };
    if needle.is_empty() {
        return best;
    }

    for rec in haystack {
        // v is the automaton state of the current match, l its length:
        // the longest suffix of the reference scanned so far that also
        // occurs in the region. A base that cannot extend the match
        // falls back along the suffix links, keeping the longest suffix
        // that still can.
        let mut v = 0usize;
        let mut l = 0u64;
        for (i, &b_raw) in rec.seq.iter().enumerate() {
            let b = b_raw.to_ascii_uppercase();
            while v != 0 && !sam.states[v].next.contains_key(&b) {
                v = sam.states[v].link;
                l = sam.states[v].len;
            }
            match sam.states[v].next.get(&b) {
                Some(&next) => {
                    v = next;
                    l += 1;
                }
                None => {
                    // Not even from the root: this base never occurs in
                    // the region. Start over from nothing.
                    v = 0;
                    l = 0;
                }
            }
            if l > best.length {
                // Every string of a state shares one endpos set, so the
                // state's recorded end position is a valid occurrence
                // end of the matched string - the example the caller
                // gets to see.
                let qry_end = sam.states[v].end_pos + 1;
                best = LongestMatch {
                    length: l,
                    ref_seqid: rec.id.clone(),
                    ref_start: (i as u64) + 2 - l,
                    ref_end: (i as u64) + 1,
                    qry_start: qry_end + 1 - l,
                    qry_end,
                };
            }
        }
    }
    best
}

/// A suffix automaton over the needle's bytes, uppercased. Standard
/// construction (Blumer et al. / Crochemore): each state is an
/// equivalence class of substrings sharing one endpos set, with the
/// suffix link pointing at the class of the longest proper suffix.
struct SuffixAutomaton {
    states: Vec<State>,
}

struct State {
    /// Length of the longest substring in this class; its suffixes down
    /// to (but excluding) the link's length all belong here too.
    len: u64,
    /// Suffix link: the class of the longest proper suffix of those.
    link: usize,
    /// Transitions, one per following base.
    next: std::collections::HashMap<u8, usize>,
    /// 0-based end position of an occurrence of this class's longest
    /// substring in the needle. Any occurrence serves the caller, which
    /// only wants one example.
    end_pos: u64,
}

impl SuffixAutomaton {
    fn new(s: &[u8]) -> SuffixAutomaton {
        let mut sam = SuffixAutomaton {
            states: vec![State {
                len: 0,
                link: 0,
                next: std::collections::HashMap::new(),
                end_pos: 0,
            }],
        };
        let mut last = 0usize;
        for (i, &b) in s.iter().enumerate() {
            let cur = sam.new_state();
            sam.states[cur].len = sam.states[last].len + 1;
            sam.states[cur].end_pos = i as u64;

            // Walk up from last, adding the transition to cur on every
            // state that lacks one. `p` as an Option is the sentinel the
            // textbook -1: None means the walk fell off the root, which
            // is exactly the case where b has never been seen before and
            // cur's suffix link is the root.
            let mut p = Some(last);
            while let Some(pp) = p {
                if let std::collections::hash_map::Entry::Vacant(e) = sam.states[pp].next.entry(b) {
                    e.insert(cur);
                    p = if pp == 0 {
                        None
                    } else {
                        Some(sam.states[pp].link)
                    };
                } else {
                    break;
                }
            }
            match p {
                None => {
                    sam.states[cur].link = 0;
                }
                Some(pp) => {
                    let q = sam.states[pp].next[&b];
                    if sam.states[pp].len + 1 == sam.states[q].len {
                        sam.states[cur].link = q;
                    } else {
                        // q's class must split: the strings that are
                        // proper suffixes of longest(pp)+b get their own
                        // state, the clone, which inherits q's outgoing
                        // transitions (any continuation of those
                        // suffixes is one of q's continuations too).
                        let clone = sam.new_state();
                        sam.states[clone].len = sam.states[pp].len + 1;
                        sam.states[clone].next = sam.states[q].next.clone();
                        sam.states[clone].link = sam.states[q].link;
                        sam.states[clone].end_pos = sam.states[q].end_pos;
                        let mut r = Some(pp);
                        while let Some(rp) = r {
                            if sam.states[rp].next.get(&b) == Some(&q) {
                                sam.states[rp].next.insert(b, clone);
                                r = if rp == 0 {
                                    None
                                } else {
                                    Some(sam.states[rp].link)
                                };
                            } else {
                                break;
                            }
                        }
                        sam.states[q].link = clone;
                        sam.states[cur].link = clone;
                    }
                }
            }
            last = cur;
        }
        sam
    }

    fn new_state(&mut self) -> usize {
        self.states.push(State {
            len: 0,
            link: 0,
            next: std::collections::HashMap::new(),
            end_pos: 0,
        });
        self.states.len() - 1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rec(id: &str, seq: &str) -> FastaRecord {
        FastaRecord {
            id: id.into(),
            desc: String::new(),
            seq: seq.as_bytes().to_vec(),
        }
    }

    #[test]
    fn finds_the_longest_common_substring() {
        let m = longest_exact_match(b"ACGTACGTAC", &[rec("chr1", "TTACGTAA")]);
        assert_eq!(m.length, 6, "TACGTA is shared");
        assert_eq!(m.ref_seqid, "chr1");
        assert_eq!((m.ref_start, m.ref_end), (2, 7), "TTACGTAA[1..7]");
        assert_eq!((m.qry_start, m.qry_end), (4, 9), "ACGTACGTAC[3..9]");
    }

    #[test]
    fn a_match_shorter_than_the_states_longest_still_reports_its_own_spot() {
        // The shared GTTTA is shorter than its automaton state's longest
        // substring (TTTTA occurs in the region), so the region
        // coordinates must come from the match length, not the state's.
        let m = longest_exact_match(b"TTTTAGGGTTTA", &[rec("c", "GTTTAX")]);
        assert_eq!(m.length, 5, "GTTTA is shared");
        assert_eq!((m.qry_start, m.qry_end), (8, 12), "TTTTAGGGTTTA[7..12]");
        assert_eq!((m.ref_start, m.ref_end), (1, 5), "GTTTAX[0..5]");
    }

    #[test]
    fn case_is_ignored() {
        let m = longest_exact_match(b"acgt", &[rec("c", "ACGT")]);
        assert_eq!(m.length, 4);
    }

    #[test]
    fn nothing_shared_reports_zero() {
        let m = longest_exact_match(b"AAAA", &[rec("c", "CCC")]);
        assert_eq!(m.length, 0);
        assert_eq!(m.ref_seqid, "");
        assert_eq!((m.ref_start, m.qry_start), (0, 0));
    }

    #[test]
    fn the_best_match_across_contigs_wins() {
        let m = longest_exact_match(
            b"GGGGATCGATC",
            &[rec("c1", "ATCGA"), rec("c2", "GGGGATCGATCCCC")],
        );
        assert_eq!(m.ref_seqid, "c2");
        assert_eq!(m.length, 11);
        assert_eq!(m.ref_start, 1);
    }

    #[test]
    fn n_is_just_another_letter() {
        // An N run in the region and an N run in the reference match -
        // honest, and the caller's job to read as noise.
        let m = longest_exact_match(b"ACNNNN", &[rec("c", "TTNNNNG")]);
        assert_eq!(m.length, 4);
    }

    /// A random-input cross-check against brute force: the automaton's
    /// whole point is correctness on inputs no hand-written case covers.
    #[test]
    fn agrees_with_brute_force_on_random_inputs() {
        fn brute(needle: &[u8], hay: &[u8]) -> u64 {
            let mut best = 0;
            for i in 0..needle.len() {
                for j in 0..hay.len() {
                    let mut l = 0;
                    while i + l < needle.len()
                        && j + l < hay.len()
                        && needle[i + l].eq_ignore_ascii_case(&hay[j + l])
                    {
                        l += 1;
                    }
                    best = best.max(l);
                }
            }
            best as u64
        }
        // xorshift, so the test needs no rand crate and no network.
        let mut seed = 0x9E3779B97F4A7C15u64;
        let mut next = move || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed
        };
        for _ in 0..200 {
            let n = (next() % 40 + 5) as usize;
            let m = (next() % 60 + 5) as usize;
            let needle: Vec<u8> = (0..n)
                .map(|_| match next() % 4 {
                    0 => b'A',
                    1 => b'C',
                    2 => b'G',
                    _ => b'T',
                })
                .collect();
            let hay: Vec<u8> = (0..m)
                .map(|_| match next() % 5 {
                    0 => b'A',
                    1 => b'C',
                    2 => b'G',
                    3 => b'T',
                    _ => b'N',
                })
                .collect();
            let want = brute(&needle, &hay);
            let hay_rec = FastaRecord {
                id: "h".into(),
                desc: String::new(),
                seq: hay,
            };
            let got = longest_exact_match(&needle, &[hay_rec]).length;
            assert_eq!(got, want, "needle {:?}", String::from_utf8_lossy(&needle));
        }
    }
}
