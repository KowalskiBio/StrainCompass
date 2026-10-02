//! Variant sets for panel genes.
//!
//! A gene name does not pin down one sequence. Listeria carries several
//! cadA genes (Tn5422, pLM80, LGI2...), around 70 % identical as protein
//! and far less as DNA, and a curated catalog may file quite different
//! ones under one name. A panel searched with the wrong one reports the
//! gene absent where it sits in full. So a panel holds every variant it
//! can find for a gene, and the search reports which one matched.
//!
//! Besides the catalog's own entries, variants come from the records the
//! panel's other genes were taken from: cadC came from the Tn5422 record,
//! so the cadA beside it there is a cadA worth searching. Genes of one
//! operon taken from different sources is exactly how the cadA of
//! Tn5422 was missed.

use crate::fasta::FastaRecord;
use crate::panel::{variant_gene, variant_id};
use std::collections::HashSet;

/// Where a panel record's sequence sits in a GenBank record.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Origin {
    pub accession: String,
    pub lo: u64,
    pub hi: u64,
}

/// The origin a panel builder note names: the catalog's "... [L28104.1:2652-2293]"
/// or the NCBI fetch's "NCBI NZ_CP0001.1:100-900". None for a reference
/// gene, a VFDB entry (a protein id) or a whole-record fetch.
pub fn record_origin(desc: &str) -> Option<Origin> {
    let from_brackets = desc
        .rfind('[')
        .and_then(|i| desc[i + 1..].split(']').next())
        .and_then(parse_range);
    from_brackets.or_else(|| {
        let mut w = desc.split_whitespace();
        (w.next() == Some("NCBI"))
            .then(|| w.next())
            .flatten()
            .and_then(parse_range)
    })
}

fn parse_range(s: &str) -> Option<Origin> {
    let (acc, range) = s.trim().rsplit_once(':')?;
    let (a, b) = range.split_once('-')?;
    let (a, b): (u64, u64) = (a.trim().parse().ok()?, b.trim().parse().ok()?);
    let ok = acc.contains('.')
        && acc
            .chars()
            .all(|c| c.is_ascii_alphanumeric() || c == '_' || c == '.');
    ok.then(|| Origin {
        accession: acc.to_string(),
        lo: a.min(b),
        hi: a.max(b),
    })
}

/// One CDS of an NCBI `fasta_cds_na` download.
#[derive(Debug, Clone, PartialEq)]
pub struct Cds {
    /// The annotated gene name, "_2"-style copy suffix removed ("cadA_2"
    /// -> "cadA"); empty when the CDS has none.
    pub gene: String,
    pub accession: String,
    pub start: u64,
    pub end: u64,
    pub seq: Vec<u8>,
}

/// Parse NCBI's `rettype=fasta_cds_na` text for record `accession`.
pub fn parse_cds_fasta(text: &str, accession: &str) -> Vec<Cds> {
    let mut out: Vec<Cds> = Vec::new();
    for line in text.lines() {
        if let Some(h) = line.strip_prefix('>') {
            let tag = |k: &str| {
                h.split(&format!("[{k}="))
                    .nth(1)
                    .and_then(|r| r.split(']').next())
                    .unwrap_or("")
                    .to_string()
            };
            let loc = tag("location");
            let nums: Vec<u64> = loc
                .split(|c: char| !c.is_ascii_digit())
                .filter_map(|n| n.parse().ok())
                .collect();
            out.push(Cds {
                gene: base_gene_name(&tag("gene")).to_string(),
                accession: accession.to_string(),
                start: nums.iter().copied().min().unwrap_or(0),
                end: nums.iter().copied().max().unwrap_or(0),
                seq: Vec::new(),
            });
        } else if let Some(c) = out.last_mut() {
            c.seq
                .extend(line.trim().bytes().map(|b| b.to_ascii_uppercase()));
        }
    }
    out
}

/// "cadA_2" -> "cadA": NCBI numbers the copies of a gene in one genome.
pub fn base_gene_name(name: &str) -> &str {
    match name.rsplit_once('_') {
        Some((g, n)) if !g.is_empty() && !n.is_empty() && n.bytes().all(|b| b.is_ascii_digit()) => {
            g
        }
        _ => name,
    }
}

const K: usize = 21;

fn kmers(s: &[u8]) -> HashSet<Vec<u8>> {
    let up = s.to_ascii_uppercase();
    let rc = crate::fasta::revcomp(&up);
    let mut out = HashSet::new();
    for v in [&up, &rc] {
        for w in v.windows(K) {
            out.insert(w.to_vec());
        }
    }
    out
}

/// Two sequences are one variant when the shorter shares nine in ten of
/// its 21-mers with the other, in either orientation: a few scattered
/// SNPs, not a different gene.
pub fn near_identical(a: &[u8], b: &[u8]) -> bool {
    let (short, long) = if a.len() <= b.len() { (a, b) } else { (b, a) };
    if short.len() < K {
        return short.eq_ignore_ascii_case(long);
    }
    let lk = kmers(long);
    let windows: Vec<&[u8]> = short.windows(K).collect();
    let shared = windows
        .iter()
        .filter(|w| lk.contains(&w.to_ascii_uppercase()))
        .count();
    shared * 10 >= windows.len() * 9
}

/// The CDS around one panel record's origin: the record it came from and
/// the genes annotated near it there.
pub struct Neighbourhood {
    /// The panel gene whose record this is.
    pub source_gene: String,
    pub origin: Origin,
    /// Short description of the record (its NCBI title).
    pub title: String,
    pub cds: Vec<Cds>,
}

/// A variant to add to the panel, with the note telling the user why.
#[derive(Debug, Clone, PartialEq)]
pub struct Addition {
    pub gene: String,
    /// The FASTA record, header note included.
    pub record: String,
    pub note: String,
}

/// Panel genes found beside another panel gene's origin, in a version the
/// panel does not yet hold: added as further variants.
pub fn partner_variants(panel: &[FastaRecord], hoods: &[Neighbourhood]) -> Vec<Addition> {
    let mut genes: Vec<String> = Vec::new();
    for r in panel {
        let g = variant_gene(&r.id).to_string();
        if !genes.contains(&g) {
            genes.push(g);
        }
    }
    let mut held: Vec<(String, Vec<u8>)> = panel
        .iter()
        .map(|r| (variant_gene(&r.id).to_string(), r.seq.clone()))
        .collect();
    let mut out: Vec<Addition> = Vec::new();
    for h in hoods {
        for c in &h.cds {
            let Some(gene) = genes.iter().find(|g| g.eq_ignore_ascii_case(&c.gene)) else {
                continue;
            };
            if c.seq.len() < 90
                || held
                    .iter()
                    .any(|(g, s)| g == gene && near_identical(s, &c.seq))
            {
                continue;
            }
            let n = held.iter().filter(|(g, _)| g == gene).count() + 1;
            let id = variant_id(gene, n);
            let mut record = format!(
                ">{id} NCBI {}:{}-{} beside {} in {}\n",
                c.accession, c.start, c.end, h.source_gene, c.accession
            );
            for chunk in c.seq.chunks(60) {
                record.push_str(&String::from_utf8_lossy(chunk));
                record.push('\n');
            }
            let first_from = panel
                .iter()
                .find(|r| r.id == *gene)
                .map(|r| match record_origin(&r.desc) {
                    Some(o) => o.accession,
                    None if r.desc.starts_with("reference") || r.desc.is_empty() => {
                        "your reference genome".into()
                    }
                    None => "another source".into(),
                })
                .unwrap_or_default();
            let note = if *gene == h.source_gene {
                format!(
                    "{gene}: also searching the {gene} annotated in {} ({}).",
                    c.accession, h.title
                )
            } else {
                format!(
                    "{gene}: {src} was taken from {acc} ({title}) but {gene} from {first_from}; genes of one operon should come from the same source, so the {gene} beside {src} in {acc} is searched as another variant.",
                    src = h.source_gene,
                    acc = c.accession,
                    title = h.title,
                )
            };
            held.push((gene.clone(), c.seq.clone()));
            out.push(Addition {
                gene: gene.clone(),
                record,
                note,
            });
        }
    }
    out
}

/// One line per gene with more than one variant, naming where each comes
/// from, so a surprising choice is visible before the run.
pub fn variant_summary(panel: &[FastaRecord]) -> Vec<String> {
    let mut genes: Vec<&str> = Vec::new();
    for r in panel {
        let g = variant_gene(&r.id);
        if !genes.contains(&g) {
            genes.push(g);
        }
    }
    genes
        .into_iter()
        .filter_map(|g| {
            let from: Vec<String> = panel
                .iter()
                .filter(|r| variant_gene(&r.id) == g)
                .map(|r| match record_origin(&r.desc) {
                    Some(o) => o.accession,
                    None if r.desc.is_empty() || r.desc.starts_with("reference") => {
                        "reference genome".into()
                    }
                    None => r.desc.split_whitespace().take(2).collect::<Vec<_>>().join(" "),
                })
                .collect();
            (from.len() > 1).then(|| {
                format!(
                    "{g}: {} sequences are searched as variants ({}); each strain's result says which one matched.",
                    from.len(),
                    from.join(", ")
                )
            })
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rec(id: &str, desc: &str, seq: &str) -> FastaRecord {
        FastaRecord {
            id: id.into(),
            desc: desc.into(),
            seq: seq.as_bytes().to_vec(),
        }
    }

    // two unrelated 120 bp genes
    const A: &str = "ATGGCTAAAGAAACTGTTTATCGTGTAGACGGTTTATCTTGTACGAATTGTGCAGCTAAATTTGAGCGAAATGTGAAAGAAATCGAAGGCGTAACTGAAGCAATCGTTAATTTTGGG";
    const B: &str = "ATGAGTAAAGCGAGTAAACAAACAACATACCGTGTAGACGGTATGTCTTGTACGAATTGTGCAGGAAAGTTTGAAAAAAACGTAAAACAACTTGCAGGAGTCCAAGATGCAAAAGTA";

    #[test]
    fn reads_origins_from_builder_notes() {
        let o = record_origin("AMRFinderPlus cadC_Lm: CadC [L28104.1:2652-2293]").unwrap();
        assert_eq!((o.accession.as_str(), o.lo, o.hi), ("L28104.1", 2293, 2652));
        let o = record_origin("NCBI NZ_CP045969.1:100-900").unwrap();
        assert_eq!(o.accession, "NZ_CP045969.1");
        assert_eq!(
            record_origin("VFDB llsG: LlsG [VFG045332(gb|AHK25017) [Listeria innocua]]"),
            None
        );
        assert_eq!(record_origin("reference lmo0200"), None);
        assert_eq!(record_origin("NCBI HF565366.1"), None);
    }

    #[test]
    fn parses_ncbi_cds_downloads() {
        let t = ">lcl|L28104.1_cds_AAA25275.1_1 [gene=cadA] [protein=ATPase] [location=complement(158..2293)] [gbkey=CDS]\nATGAAA\nTAA\n\
                 >lcl|AP022822.1_cds_X_2 [gene=cadA_2] [location=2608314..2610431]\nATG\n\
                 >lcl|X_3 [protein=hypothetical] [location=1..9]\nATG\n";
        let c = parse_cds_fasta(t, "L28104.1");
        assert_eq!(c.len(), 3);
        assert_eq!(
            (c[0].gene.as_str(), c[0].start, c[0].end),
            ("cadA", 158, 2293)
        );
        assert_eq!(c[0].seq, b"ATGAAATAA");
        assert_eq!(c[1].gene, "cadA");
        assert_eq!(c[2].gene, "");
    }

    #[test]
    fn near_identical_tolerates_snps_not_other_genes() {
        let mut snp = A.as_bytes().to_vec();
        snp[60] = b'T';
        assert!(near_identical(A.as_bytes(), &snp));
        assert!(near_identical(
            A.as_bytes(),
            &crate::fasta::revcomp(A.as_bytes())
        ));
        assert!(!near_identical(A.as_bytes(), B.as_bytes()));
    }

    #[test]
    fn adds_the_operon_partner_from_the_other_genes_record() {
        // cadA from the catalog's Enterococcus entry, cadC from Tn5422
        let panel = vec![
            rec(
                "cadA",
                "AMRFinderPlus cadA_Lm: CadA [AP022822.1:2608314-2610431]",
                B,
            ),
            rec(
                "cadC",
                "AMRFinderPlus cadC_Lm: CadC [L28104.1:2652-2293]",
                "ATGCCC",
            ),
        ];
        let hood = Neighbourhood {
            source_gene: "cadC".into(),
            origin: record_origin(&panel[1].desc).unwrap(),
            title: "Listeria monocytogenes transposon Tn5422".into(),
            cds: vec![
                Cds {
                    gene: "cadA".into(),
                    accession: "L28104.1".into(),
                    start: 158,
                    end: 2293,
                    seq: A.as_bytes().to_vec(),
                },
                Cds {
                    gene: "tnpR".into(),
                    accession: "L28104.1".into(),
                    start: 2932,
                    end: 3486,
                    seq: A.as_bytes().to_vec(),
                },
            ],
        };
        let add = partner_variants(&panel, std::slice::from_ref(&hood));
        assert_eq!(add.len(), 1);
        assert_eq!(add[0].gene, "cadA");
        assert!(add[0]
            .record
            .starts_with(">cadA__v2 NCBI L28104.1:158-2293 beside cadC"));
        assert!(add[0].note.contains("AP022822.1") && add[0].note.contains("same source"));
        // a second pass over the same record adds nothing new
        let mut panel2 = panel.clone();
        panel2.push(rec("cadA__v2", "NCBI L28104.1:158-2293", A));
        assert!(partner_variants(&panel2, &[hood]).is_empty());
    }

    #[test]
    fn summarises_genes_with_several_variants() {
        let panel = vec![
            rec("cadA", "AMRFinderPlus cadA_Lm: x [AP022822.1:1-9]", A),
            rec("cadA__v2", "NCBI L28104.1:158-2293 beside cadC", B),
            rec("inlA", "reference lmo0433", A),
        ];
        let s = variant_summary(&panel);
        assert_eq!(s.len(), 1);
        assert!(s[0].starts_with("cadA: 2 sequences") && s[0].contains("AP022822.1, L28104.1"));
    }
}
