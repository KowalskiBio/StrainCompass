//! The local reference library (see `straincompass_engine::library`) as
//! the app uses it: where libraries are installed, which one a project
//! uses (the one of its genus), and the panel lookups built on it.
//!
//! Installed under STRAINCOMPASS_LIBRARY_DIR, else `library/` beside the
//! data directory - never inside it, since every deploy backup copies the
//! whole data directory.

use crate::error::ApiResult;
use crate::state::SharedState;
use axum::extract::State;
use axum::Json;
use serde::Serialize;
use std::path::{Path, PathBuf};
use straincompass_engine::library::{Library, LibraryVariant};
use straincompass_engine::panel::variant_id;

pub fn root(state: &SharedState) -> PathBuf {
    std::env::var_os("STRAINCOMPASS_LIBRARY_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            state
                .data_dir
                .parent()
                .map(|p| p.join("library"))
                .unwrap_or_else(|| PathBuf::from("library"))
        })
}

/// The library of `genus`, if one is installed and readable. A broken one
/// is logged and treated as absent, so the app falls back to NCBI rather
/// than failing; GET /library shows the error.
pub fn for_genus(state: &SharedState, genus: &str) -> Option<Library> {
    match Library::for_genus(&root(state), genus)? {
        Ok(lib) => Some(lib),
        Err(e) => {
            tracing::warn!("reference library for {genus}: {e}");
            None
        }
    }
}

/// Every usable installed library, one per genus.
pub fn installed(state: &SharedState) -> Vec<Library> {
    let Ok(dirs) = std::fs::read_dir(root(state)) else {
        return Vec::new();
    };
    let mut out: Vec<Library> = dirs
        .filter_map(|d| d.ok())
        .filter_map(|d| {
            let genus = d.file_name().to_string_lossy().to_string();
            for_genus(state, &genus)
        })
        .collect();
    out.sort_by(|a, b| a.manifest.genus.cmp(&b.manifest.genus));
    out
}

/// A genome stored in any installed library, by assembly accession: the
/// reference step runs before the project has a genome to read a genus
/// from, so every library is asked.
pub fn find_reference(
    state: &SharedState,
    accession: &str,
) -> Option<(String, straincompass_engine::library::StoredReference)> {
    installed(state).into_iter().find_map(|lib| {
        lib.reference(accession)
            .ok()
            .flatten()
            .map(|r| (lib.title(), r))
    })
}

pub fn for_project(state: &SharedState, project_id: i64) -> Option<Library> {
    let genus = crate::routes::origin::project_genus(state, project_id)?;
    for_genus(state, &genus)
}

#[derive(Serialize)]
pub struct LibraryStatus {
    pub genus: String,
    pub version: String,
    pub built: String,
    pub counts: std::collections::BTreeMap<String, u64>,
    /// Why an installed library cannot be used, if it cannot.
    pub error: Option<String>,
}

#[derive(Serialize)]
pub struct LibrariesDto {
    pub root: String,
    pub libraries: Vec<LibraryStatus>,
}

/// GET /library : the installed libraries, one per genus.
pub async fn list(State(state): State<SharedState>) -> ApiResult<Json<LibrariesDto>> {
    let root = root(&state);
    let mut libraries = Vec::new();
    if let Ok(dirs) = std::fs::read_dir(&root) {
        let mut dirs: Vec<PathBuf> = dirs.filter_map(|d| d.ok().map(|d| d.path())).collect();
        dirs.sort();
        for d in dirs {
            let current = d.join("current");
            if !current.join("manifest.json").is_file() {
                continue;
            }
            let genus = d
                .file_name()
                .unwrap_or_default()
                .to_string_lossy()
                .to_string();
            libraries.push(match Library::open(&current) {
                Ok(lib) => LibraryStatus {
                    genus: lib.manifest.genus,
                    version: lib.manifest.version,
                    built: lib.manifest.built,
                    counts: lib.manifest.counts,
                    error: None,
                },
                Err(e) => LibraryStatus {
                    genus,
                    version: String::new(),
                    built: String::new(),
                    counts: Default::default(),
                    error: Some(e.to_string()),
                },
            });
        }
    }
    Ok(Json(LibrariesDto {
        root: root.display().to_string(),
        libraries,
    }))
}

/// What the library contributed to a panel.
#[derive(Default)]
pub struct PanelLookup {
    /// Panel records to append.
    pub records: String,
    /// One line per gene taken from the library, for the user.
    pub notes: Vec<String>,
    /// Further variants of reference genes, in words.
    pub hints: Vec<String>,
    /// The wanted lines the library does not know.
    pub unresolved: Vec<String>,
}

/// Resolve the panel's genes against the library: every variant of each
/// wanted name not in the reference, and for each reference gene the
/// variants other than its own. Lines pinning a GenBank record are the
/// user's explicit choice and are left alone.
///
/// `blastn` is needed only to tell a reference gene's own group from its
/// other variants; without it reference genes get no extra variants.
/// `blastx` keeps the variants to genes related to the one meant (a name
/// can belong to unrelated genes: prfA is the PrfA regulator and release
/// factor 1); without it every gene of the name is searched.
#[allow(clippy::too_many_arguments)]
pub fn panel_lookup(
    lib: &Library,
    blastn: Option<&Path>,
    blastx: Option<&Path>,
    wanted: &[String],
    found: &[String],
    ref_panel: &str,
    ref_fasta: Option<&Path>,
    work: &Path,
) -> PanelLookup {
    let mut out = PanelLookup::default();
    let title = lib.title();
    for line in wanted {
        if let Some(spec) = crate::routes::ncbi::parse_accession_spec(line) {
            // a pinned record: the library's copy of it, else NCBI's
            match pinned(lib, &spec) {
                Some((record, note)) => {
                    out.records.push_str(&record);
                    out.notes.push(note);
                }
                None => out.unresolved.push(line.clone()),
            }
            continue;
        }
        let id = line
            .split_whitespace()
            .next()
            .unwrap_or("")
            .trim_matches(|c| c == '(' || c == ')' || c == '[' || c == ']')
            .to_string();
        // "pva (lmo0446)": the locus tag first, then the symbol
        let variants = lookup_tokens(line)
            .iter()
            .find_map(|t| lib.variants(t).ok().filter(|v| !v.is_empty()));
        let Some(variants) = variants.filter(|_| !id.is_empty()) else {
            out.unresolved.push(line.clone());
            continue;
        };
        // the first variant is the one a curated entry names, or else the
        // most common: the others must be related to it
        let anchor = variants[0].seq.clone();
        let products = vec![variants[0].product.clone()];
        let (variants, unrelated) =
            split_family(lib, blastx, &anchor, &products, variants, Some(0), work);
        if !unrelated.is_empty() {
            out.hints.push(unrelated_note(&id, &unrelated));
        }
        for (i, v) in variants.iter().enumerate() {
            out.records
                .push_str(&record(&variant_id(&id, i + 1), v, &title));
        }
        out.notes.push(format!(
            "{id} \u{2190} {}",
            variants.iter().map(summary).collect::<Vec<_>>().join("; ")
        ));
    }

    // reference genes: the library's other variants of each, beside the
    // reference copy (which stays the first)
    let Some(blastn) = blastn else {
        return out;
    };
    let Ok(recs) = straincompass_engine::fasta::parse_fasta_str(ref_panel) else {
        return out;
    };
    let own: Vec<(String, Vec<u8>)> = found
        .iter()
        .filter_map(|name| {
            recs.iter()
                .find(|r| &r.id == name)
                .map(|r| (name.clone(), r.seq.clone()))
        })
        .collect();
    let own_groups = match lib.groups_of(blastn, &own, work) {
        Ok(g) => g,
        Err(e) => {
            tracing::warn!("reference library: {e}");
            return out;
        }
    };
    for (name, seq) in &own {
        let mine = own_groups.get(name).cloned().unwrap_or_default();
        let others: Vec<LibraryVariant> = lib
            .variants(name)
            .unwrap_or_default()
            .into_iter()
            .filter(|v| !mine.contains(&v.group) && !v.seq.eq_ignore_ascii_case(seq))
            .collect();
        // what the reference copy is, by its own groups' products
        let products: Vec<String> = mine
            .iter()
            .filter_map(|g| lib.group_brief(*g).ok().flatten().map(|b| b.2))
            .collect();
        let (others, unrelated) = split_family(lib, blastx, seq, &products, others, None, work);
        if !unrelated.is_empty() {
            out.hints.push(unrelated_note(name, &unrelated));
        }
        // a variant the reference itself carries elsewhere is a paralog
        // of the reference gene (EGD-e's lmo1490 beside its aroE, lmo0490)
        let seqs: Vec<Vec<u8>> = others.iter().map(|v| v.seq.clone()).collect();
        let carried = match ref_fasta.map(|f| lib.carried_by(blastn, &seqs, f, work)) {
            Some(Ok(c)) => c,
            Some(Err(e)) => {
                tracing::warn!("reference library: {e}");
                vec![false; others.len()]
            }
            None => vec![false; others.len()],
        };
        let (others, paralogs): (Vec<_>, Vec<_>) =
            others.into_iter().zip(carried).partition(|(_, c)| !c);
        let others: Vec<LibraryVariant> = others.into_iter().map(|(v, _)| v).collect();
        if !paralogs.is_empty() {
            out.hints.push(format!(
                "{name}: your reference also carries {} other gene{} the library files under {name} ({}): {}, not other versions of it, so {} not searched as {name}.",
                paralogs.len(),
                if paralogs.len() == 1 { "" } else { "s" },
                paralogs
                    .iter()
                    .map(|(v, _)| if v.product.is_empty() { v.label.clone() } else { v.product.clone() })
                    .collect::<Vec<_>>()
                    .join("; "),
                if paralogs.len() == 1 { "a paralog" } else { "paralogs" },
                if paralogs.len() == 1 { "it is" } else { "they are" },
            ));
        }
        if others.is_empty() {
            continue;
        }
        for (i, v) in others.iter().enumerate() {
            // numbered after the reference copy; renumbered with the rest
            out.records
                .push_str(&record(&variant_id(name, i + 2), v, &title));
        }
        out.hints.push(format!(
            "{name}: taken from your reference genome. The {title} holds {} other variant{} of {name} ({}); {} searched too, and the results say which one matched.",
            others.len(),
            if others.len() == 1 { "" } else { "s" },
            others.iter().map(summary).collect::<Vec<_>>().join("; "),
            if others.len() == 1 { "it is" } else { "they are" },
        ));
    }
    out
}

/// Largest whole record taken as a panel gene (a transposon or cassette),
/// as for NCBI records: whole chromosomes are not genes.
const MAX_WHOLE_RECORD: usize = 25_000;

/// A gene pinned to a GenBank record ("qacH (HF565366.1)", "emrC
/// (CP038643.1:1496-1882 rev)"), read from the library's copy of the
/// record: (panel record, note). None when the library lacks the record,
/// or its annotation does not name the gene - NCBI's GenBank original
/// may, so the line is left to it.
fn pinned(lib: &Library, spec: &crate::routes::ncbi::AccessionSpec) -> Option<(String, String)> {
    use straincompass_engine::library::Replicon;
    let rep = lib.replicon(&spec.accession).ok()??;
    let seq = lib.replicon_seq(&rep).ok()??;
    let name = spec.name.clone().unwrap_or_else(|| spec.accession.clone());
    let (lo, hi, minus, product) = match (spec.coords, &spec.name) {
        (Some((s, e, rev)), _) => (s.min(e).max(1) as u64, s.max(e) as u64, rev, String::new()),
        (None, Some(n)) => {
            let base =
                |s: &str| straincompass_engine::panel_variants::base_gene_name(s).to_string();
            let genes = lib.replicon_genes(&rep).ok()?;
            let g = genes.iter().find(|g| {
                [
                    base(&g.symbol),
                    g.locus_tag.clone(),
                    g.old_locus_tag.clone(),
                ]
                .iter()
                .any(|t| !t.is_empty() && t.eq_ignore_ascii_case(n))
            })?;
            (g.start, g.end, g.strand < 0, g.product.clone())
        }
        (None, None) if seq.len() <= MAX_WHOLE_RECORD => {
            (1, seq.len() as u64, false, String::new())
        }
        (None, None) => return None,
    };
    let gene = Replicon::slice(&seq, lo, hi, minus)?;
    let place = format!(
        "{}:{lo}-{hi}{}",
        rep.accession,
        if minus { " rev" } else { "" }
    );
    let what = if product.is_empty() {
        rep.title()
    } else {
        product
    };
    let mut record = format!(
        ">{name} Library {}: {what} [{}:{}-{}]\n",
        rep.accession,
        rep.accession,
        if minus { hi } else { lo },
        if minus { lo } else { hi },
    );
    for chunk in gene.chunks(60) {
        record.push_str(&String::from_utf8_lossy(chunk));
        record.push('\n');
    }
    Some((
        record,
        format!(
            "{name} \u{2190} {place}, the {} copy of {} ({})",
            lib.title(),
            spec.accession,
            rep.title()
        ),
    ))
}

/// `variants` split into those related to `anchor` (a gene's DNA) and
/// those sharing only the name. `keep` is always related (the anchor's
/// own variant). Without blastx, or when the search fails, all count as
/// related, as before.
fn split_family(
    lib: &Library,
    blastx: Option<&Path>,
    anchor: &[u8],
    products: &[String],
    variants: Vec<LibraryVariant>,
    keep: Option<usize>,
    work: &Path,
) -> (Vec<LibraryVariant>, Vec<LibraryVariant>) {
    // the same product is the same gene even past what blastx links (the
    // ActA of L. ivanovii and of L. monocytogenes)
    let same_product = |v: &LibraryVariant| {
        let p = v.product.trim().to_ascii_lowercase();
        !p.is_empty()
            && p != "hypothetical protein"
            && products.iter().any(|q| q.trim().eq_ignore_ascii_case(&p))
    };
    let family = match blastx.map(|b| lib.family(b, anchor, work)) {
        Some(Ok(f)) => f,
        Some(Err(e)) => {
            tracing::warn!("reference library: {e}");
            return (variants, Vec::new());
        }
        None => return (variants, Vec::new()),
    };
    let (mut related, mut unrelated) = (Vec::new(), Vec::new());
    for (i, v) in variants.into_iter().enumerate() {
        if Some(i) == keep || family.contains(&v.group) || same_product(&v) {
            related.push(v);
        } else {
            unrelated.push(v);
        }
    }
    (related, unrelated)
}

/// "prfA: the Listeria library also files 1 unrelated gene under prfA
/// (peptide chain release factor 1 on chromosome (...), 832 genomes); it
/// is not searched as prfA. To search it, add it under a name of its own,
/// e.g. prfA_2 (NZ_CP...:1-1077 rev)."
fn unrelated_note(name: &str, unrelated: &[LibraryVariant]) -> String {
    let first = &unrelated[0];
    format!(
        "{name}: the reference library also files {} unrelated gene{} under {name} ({}); {} not searched as {name}. To search one, add it under a name of its own, e.g. {name}_2 ({}{}).",
        unrelated.len(),
        if unrelated.len() == 1 { "" } else { "s" },
        unrelated
            .iter()
            .map(|v| format!(
                "{} on {}, {} genome{}",
                if v.product.is_empty() { &v.label } else { &v.product },
                v.place,
                v.n_genomes,
                if v.n_genomes == 1 { "" } else { "s" }
            ))
            .collect::<Vec<_>>()
            .join("; "),
        if unrelated.len() == 1 { "it is" } else { "they are" },
        first.locus,
        if first.minus { " rev" } else { "" },
    )
}

/// The names to try for one line of a gene list, the parenthesized locus
/// tag first: "pva (lmo0446)" -> ["lmo0446", "pva"].
fn lookup_tokens(line: &str) -> Vec<String> {
    let mut tokens: Vec<String> = Vec::new();
    for part in line.split([',', '\t', ';']) {
        for tok in part.split_whitespace() {
            let bare = tok
                .trim_start_matches(['(', '['])
                .trim_end_matches([')', ']']);
            if bare.is_empty() || bare.starts_with('#') {
                continue;
            }
            if tok.starts_with('(') || tok.starts_with('[') {
                tokens.insert(0, bare.to_string());
            } else {
                tokens.push(bare.to_string());
            }
        }
    }
    tokens
}

/// A library variant as a panel record. The bracket holds more than an
/// origin range on purpose: the NCBI partner search must not refetch
/// what the library already settled.
fn record(id: &str, v: &LibraryVariant, title: &str) -> String {
    let mut rec = format!(
        ">{id} Library g{} {}: {} [{title}; {}, {}; {}; in {} genome{}]\n",
        v.group,
        v.label,
        v.product,
        v.place,
        v.locus,
        v.evidence,
        v.n_genomes,
        if v.n_genomes == 1 { "" } else { "s" },
    );
    for chunk in v.seq.chunks(60) {
        rec.push_str(&String::from_utf8_lossy(chunk));
        rec.push('\n');
    }
    rec
}

/// "plasmid pLmN1546 (Listeria monocytogenes X), 6 genomes (GenBank L28104.1 cadA)"
fn summary(v: &LibraryVariant) -> String {
    let why = match v.matched_by.as_str() {
        "homolog" => format!(", unnamed, {:.0} % like a named {}", v.identity, v.label),
        "catalog" => format!(", {}", v.evidence.trim_start_matches("matches ")),
        _ => String::new(),
    };
    format!(
        "{} on {}, {} genome{}{why}",
        v.label,
        v.place,
        v.n_genomes,
        if v.n_genomes == 1 { "" } else { "s" }
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tries_the_locus_tag_first() {
        assert_eq!(lookup_tokens("pva (lmo0446)"), ["lmo0446", "pva"]);
        assert_eq!(lookup_tokens("cadA"), ["cadA"]);
        assert_eq!(lookup_tokens("inlA\tlmo0433"), ["inlA", "lmo0433"]);
    }

    #[test]
    fn library_records_are_not_taken_for_ncbi_origins() {
        let v = LibraryVariant {
            group: 11,
            label: "cadA".into(),
            product: "heavy metal translocating P-type ATPase".into(),
            matched_by: "catalog".into(),
            identity: 100.0,
            evidence: "matches GenBank L28104.1 cadA (100 % protein identity)".into(),
            seq: b"ATGCCC".to_vec(),
            locus: "NZ_CP2.1:47228-49363".into(),
            place: "plasmid pLmN1546 (Listeria monocytogenes X)".into(),
            n_genomes: 6,
            n_plasmid: 4,
            n_chromosome: 2,
            minus: true,
        };
        let rec = record("cadA__v2", &v, "Listeria library t");
        let recs = straincompass_engine::fasta::parse_fasta_str(&rec).unwrap();
        assert_eq!(recs[0].id, "cadA__v2");
        assert_eq!(
            straincompass_engine::panel_variants::record_origin(&recs[0].desc),
            None
        );
        assert_eq!(
            straincompass_engine::element::panel_gene_source("cadA", &recs[0].desc, &[]),
            "reference library (g11)"
        );
        assert_eq!(
            summary(&v),
            "cadA on plasmid pLmN1546 (Listeria monocytogenes X), 6 genomes, GenBank L28104.1 cadA (100 % protein identity)"
        );
    }
}
