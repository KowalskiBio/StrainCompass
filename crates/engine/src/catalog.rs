//! Gene sequences by name from curated catalogs: the AMRFinderPlus
//! reference CDS set (`AMR_CDS.fa`: resistance, biocide, metal, stress)
//! and the VFDB core set (verified virulence genes). Used when a panel
//! gene is not in the reference genome, before any free-text NCBI search.
//!
//! Gene names are not unique - "emrC" is an E. coli outer membrane
//! protein and a Listeria plasmid efflux pump - and both catalogs carry
//! organism-specific entries for exactly that reason (AMRFinderPlus
//! appends an organism code: emrC_Lis, cadC_Lm, cadC_Sa). The lookup
//! prefers the entry that matches the project's organism and reports
//! which entry and product it took, so a surprising choice is visible.

use std::path::Path;

/// One catalog gene picked for a requested name.
#[derive(Debug, Clone, PartialEq)]
pub struct CatalogGene {
    /// "AMRFinderPlus" or "VFDB".
    pub source: String,
    /// The catalog's own name for it (emrC_Lis, SMR_efflux_bcrB, llsG).
    pub symbol: String,
    pub product: String,
    /// Where the sequence comes from (accession:range, or the VFDB id and
    /// its source organism).
    pub origin: String,
    pub seq: Vec<u8>,
}

/// A parsed catalog: (symbol, product, origin, organism hint, sequence).
pub struct Catalog {
    source: &'static str,
    entries: Vec<Entry>,
}

struct Entry {
    symbol: String,
    product: String,
    origin: String,
    /// VFDB: the source organism; AMRFinderPlus: empty (the organism code
    /// is in the symbol suffix instead).
    organism: String,
    seq: Vec<u8>,
}

fn records(bytes: &[u8]) -> Vec<(String, Vec<u8>)> {
    let text = String::from_utf8_lossy(bytes);
    let mut out: Vec<(String, Vec<u8>)> = Vec::new();
    for line in text.lines() {
        if let Some(h) = line.strip_prefix('>') {
            out.push((h.trim().to_string(), Vec::new()));
        } else if let Some(last) = out.last_mut() {
            last.1.extend(line.trim().bytes().map(|b| b.to_ascii_uppercase()));
        }
    }
    out
}

impl Catalog {
    /// AMRFinderPlus `AMR_CDS.fa`: `>prot|nuc|1|1|symbol|symbol|product_words acc:range`.
    pub fn amrfinder(bytes: &[u8]) -> Catalog {
        let entries = records(bytes)
            .into_iter()
            .filter_map(|(h, seq)| {
                let (fields, origin) = h.split_once(' ').unwrap_or((h.as_str(), ""));
                let f: Vec<&str> = fields.split('|').collect();
                let symbol = f.get(4)?.to_string();
                Some(Entry {
                    product: f.get(6).unwrap_or(&"").replace('_', " "),
                    origin: origin.trim().to_string(),
                    organism: String::new(),
                    symbol,
                    seq,
                })
            })
            .collect();
        Catalog {
            source: "AMRFinderPlus",
            entries,
        }
    }

    /// VFDB: `>VFGid(gb|acc) (gene/alias) product [factor - category] [organism]`.
    pub fn vfdb(bytes: &[u8]) -> Catalog {
        let entries = records(bytes)
            .into_iter()
            .filter_map(|(h, seq)| {
                let (id, gene, product, _, _, organism) = crate::screen::parse_vfdb_header(&h);
                if gene.is_empty() {
                    return None;
                }
                Some(Entry {
                    symbol: gene,
                    product,
                    origin: format!("{id} [{organism}]"),
                    organism,
                    seq,
                })
            })
            .collect();
        Catalog {
            source: "VFDB",
            entries,
        }
    }

    pub fn from_file(path: &Path, vfdb: bool) -> Option<Catalog> {
        let bytes = std::fs::read(path).ok()?;
        Some(if vfdb {
            Catalog::vfdb(&bytes)
        } else {
            Catalog::amrfinder(&bytes)
        })
    }

    /// The best entry for `name`, with its score (higher is better), or
    /// None when the catalog has no gene of that name.
    fn best(&self, name: &str, organism: Option<&str>) -> Option<(u8, &Entry)> {
        let mut best: Option<(u8, &Entry)> = None;
        for e in &self.entries {
            let Some(score) = self.score(e, name, organism) else {
                continue;
            };
            if e.seq.len() < 50 {
                continue;
            }
            if best.is_none_or(|(s, _)| score > s) {
                best = Some((score, e));
            }
        }
        best
    }

    fn score(&self, e: &Entry, name: &str, organism: Option<&str>) -> Option<u8> {
        let (genus, species) = genus_species(organism.unwrap_or(""));
        if self.source == "VFDB" {
            // "hbp1/svpA": every alias counts
            let named = e.symbol.split('/').any(|s| s == name)
                || e.symbol.split('/').any(|s| s.eq_ignore_ascii_case(name));
            if !named {
                return None;
            }
            let mut words = e.organism.split_whitespace();
            let (g, s) = (words.next().unwrap_or(""), words.next().unwrap_or(""));
            return Some(if !genus.is_empty() && g == genus && s == species {
                4
            } else if !genus.is_empty() && g == genus {
                3
            } else {
                1
            });
        }
        // AMRFinderPlus: exact, organism-suffixed (emrC_Lis) or
        // family-prefixed (SMR_efflux_bcrB)
        if e.symbol == name {
            return Some(2);
        }
        if let Some(suffix) = e.symbol.strip_prefix(&format!("{name}_")) {
            return Some(if organism_code_matches(suffix, &genus, &species) {
                4
            } else {
                1
            });
        }
        if e.symbol.ends_with(&format!("_{name}")) {
            return Some(2);
        }
        None
    }
}

fn genus_species(organism: &str) -> (String, String) {
    let mut w = organism.split_whitespace();
    (
        w.next().unwrap_or("").to_string(),
        w.next().unwrap_or("").to_string(),
    )
}

/// AMRFinderPlus organism codes: the genus' first three letters ("Lis"
/// for Listeria) or genus + species initials ("Lm", "Sa").
fn organism_code_matches(code: &str, genus: &str, species: &str) -> bool {
    if genus.is_empty() {
        return false;
    }
    let code = code.to_ascii_lowercase();
    let g = genus.to_ascii_lowercase();
    let s = species.to_ascii_lowercase();
    code == g.chars().take(3).collect::<String>()
        || (code.len() == 2
            && code.starts_with(g.chars().next().unwrap_or(' '))
            && s.starts_with(&code[1..]))
}

/// Look `name` up in all catalogs; the best-scoring entry wins, the
/// earlier catalog on a tie.
pub fn lookup(name: &str, organism: Option<&str>, catalogs: &[Catalog]) -> Option<CatalogGene> {
    let name = name.trim();
    if name.is_empty() {
        return None;
    }
    let mut best: Option<(u8, &str, &Entry)> = None;
    for c in catalogs {
        if let Some((score, e)) = c.best(name, organism) {
            if best.is_none_or(|(s, _, _)| score > s) {
                best = Some((score, c.source, e));
            }
        }
    }
    best.map(|(_, source, e)| CatalogGene {
        source: source.to_string(),
        symbol: e.symbol.clone(),
        product: e.product.clone(),
        origin: e.origin.clone(),
        seq: e.seq.clone(),
    })
}

/// An organism-specific catalog variant of a gene the reference already
/// has under the same name, when its sequence differs: cadA from an EGD-e
/// reference vs AMRFinderPlus cadA_Lm (the Tn5422 cadA, ~70 % identical).
/// The reference copy stays the answer; this lets the user see the other.
pub fn organism_variant(
    name: &str,
    organism: Option<&str>,
    reference_seq: &[u8],
    catalogs: &[Catalog],
) -> Option<CatalogGene> {
    let g = lookup(name, organism, catalogs)?;
    let specific = g.source == "AMRFinderPlus" && g.symbol.starts_with(&format!("{name}_"));
    let same = g.seq.eq_ignore_ascii_case(reference_seq);
    (specific && !same).then_some(g)
}

#[cfg(test)]
mod tests {
    use super::*;

    const AMR: &str = ">CAQ45062.1|AM743169.1|1|1|emrC|emrC|multidrug_efflux_transporter_outer_membrane_subunit_EmrC AM743169.1:1586491-1587984\n\
        ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAAAA\n\
        >EAC4468893.1|AAAIJS010000026.1|1|1|emrC_Lis|emrC_Lis|multidrug_efflux_transporter_outer_membrane_subunit_EmrC AAAIJS010000026.1:2413-2027\n\
        TTTTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n\
        >AAA25276.1|L28104.1|1|1|cadC_Lm|cadC_Lm|Cd(II)-sensing_repressor_CadC L28104.1:2652-2293\n\
        GGGGACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n\
        >AAB59153.1|J04551.1|1|1|cadC_Sa|cadC_Sa|Cd(II)/Pb(II)/Zn(II)-sensing_repressor_CadC J04551.1:703-1071\n\
        CCCCACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n\
        >WP_003725293.1|NG_055638.1|1|1|SMR_efflux_bcrB|SMR_efflux_bcrB|quaternary_ammonium_compound_efflux_SMR_transporter_BcrB NG_055638.1:101-418\n\
        AAAAACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";

    const VF: &str = ">VFG045332(gb|AHK25017) (llsG) ABC transporter ATP-binding protein LlsG [LLS (VF0410) - Exotoxin (VFC0235)] [Listeria innocua SLCC6294]\n\
        ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n\
        >VFG000001(gb|X) (hly) listeriolysin O [LLO (VF0064) - Exotoxin (VFC0235)] [Listeria ivanovii]\n\
        AAAAACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n\
        >VFG000002(gb|Y) (hly) listeriolysin O [LLO (VF0064) - Exotoxin (VFC0235)] [Listeria monocytogenes EGD-e]\n\
        CCCCACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";

    fn cats() -> Vec<Catalog> {
        vec![Catalog::amrfinder(AMR.as_bytes()), Catalog::vfdb(VF.as_bytes())]
    }

    #[test]
    fn prefers_the_entry_for_the_projects_organism() {
        let c = cats();
        let lm = Some("Listeria monocytogenes");
        assert_eq!(lookup("emrC", lm, &c).unwrap().symbol, "emrC_Lis");
        assert_eq!(lookup("emrC", Some("Escherichia coli"), &c).unwrap().symbol, "emrC");
        assert_eq!(lookup("emrC", None, &c).unwrap().symbol, "emrC");
        assert_eq!(lookup("cadC", lm, &c).unwrap().symbol, "cadC_Lm");
        assert_eq!(lookup("cadC", Some("Staphylococcus aureus"), &c).unwrap().symbol, "cadC_Sa");
        let h = lookup("hly", lm, &c).unwrap();
        assert_eq!(h.source, "VFDB");
        assert!(h.origin.contains("EGD-e"), "{}", h.origin);
    }

    #[test]
    fn finds_family_prefixed_and_vfdb_names() {
        let c = cats();
        let b = lookup("bcrB", Some("Listeria monocytogenes"), &c).unwrap();
        assert_eq!((b.symbol.as_str(), b.origin.as_str()), ("SMR_efflux_bcrB", "NG_055638.1:101-418"));
        assert!(b.product.contains("BcrB"));
        let l = lookup("llsG", Some("Listeria monocytogenes"), &c).unwrap();
        assert_eq!(l.source, "VFDB");
        assert_eq!(l.seq.len(), 60);
        assert!(lookup("lmo0444", None, &c).is_none());
        assert!(lookup("", None, &c).is_none());
    }

    #[test]
    fn flags_a_different_organism_variant_of_a_reference_gene() {
        let c = cats();
        let lm = Some("Listeria monocytogenes");
        let v = organism_variant("cadC", lm, b"ACGT", &c).unwrap();
        assert_eq!(v.symbol, "cadC_Lm");
        // the reference already carries exactly that sequence: nothing to say
        assert!(organism_variant("cadC", lm, &v.seq, &c).is_none());
        // a plain catalog name is not an organism variant
        assert!(organism_variant("hly", lm, b"ACGT", &c).is_none());
    }

    #[test]
    fn reads_organism_codes() {
        assert!(organism_code_matches("Lis", "Listeria", "monocytogenes"));
        assert!(organism_code_matches("Lm", "Listeria", "monocytogenes"));
        assert!(organism_code_matches("Sa", "Staphylococcus", "aureus"));
        assert!(!organism_code_matches("Sa", "Listeria", "monocytogenes"));
        assert!(!organism_code_matches("Lm", "", ""));
    }
}
