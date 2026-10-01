//! Resistance and virulence screen of one query genome, independent of
//! the gene panel: AMRFinderPlus (NCBI's curated antibiotic resistance,
//! disinfectant, metal and stress genes, matched by sequence, so a name
//! shared by unrelated genes cannot mislead it) and a BLAST of the VFDB
//! core set A (experimentally verified virulence genes).
//!
//! Both resources are optional on a server: a missing one is reported as
//! such and never fails the run.

use crate::tools::{find_optional, ToolPaths};
use crate::{friendly, EngineError, Result};
use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::OnceLock;
use straincompass_types::{ScreenHit, ScreenStatus};

/// abricate's defaults for a VFDB call.
pub const VFDB_MIN_ID: f64 = 80.0;
pub const VFDB_MIN_COV: f64 = 80.0;

/// Where the screen's resources are.
#[derive(Debug, Clone, Default)]
pub struct ScreenConfig {
    pub amrfinder: Option<PathBuf>,
    /// AMRFinderPlus database directory (`-d`); its own default when unset.
    pub amr_db: Option<PathBuf>,
    /// VFDB core set A, nucleotide FASTA.
    pub vfdb: Option<PathBuf>,
    pub threads: usize,
}

impl ScreenConfig {
    /// AMRFinderPlus from the tool directories or PATH, its database from
    /// STRAINCOMPASS_AMRFINDER_DB, the VFDB file from STRAINCOMPASS_VFDB or
    /// `vfdb_default` when that file exists.
    pub fn discover(vfdb_default: Option<PathBuf>) -> ScreenConfig {
        let env_path = |k: &str| std::env::var_os(k).map(PathBuf::from);
        ScreenConfig {
            amrfinder: find_optional("amrfinder"),
            amr_db: env_path("STRAINCOMPASS_AMRFINDER_DB").filter(|p| p.is_dir()),
            vfdb: env_path("STRAINCOMPASS_VFDB")
                .or(vfdb_default)
                .filter(|p| p.is_file()),
            threads: 2,
        }
    }
}

/// Screen one query. Never errors on a missing resource: the status says
/// what ran.
pub fn screen_query(
    tools: &ToolPaths,
    cfg: &ScreenConfig,
    qry_fasta: &Path,
    organism: Option<&str>,
    work: &Path,
) -> Result<(Vec<ScreenHit>, ScreenStatus)> {
    std::fs::create_dir_all(work)?;
    let mut hits = Vec::new();
    // each half stands on its own: one failing keeps the other's genes
    let mut missing: Vec<String> = Vec::new();
    match &cfg.amrfinder {
        Some(bin) => match run_amrfinder(bin, cfg, qry_fasta, organism, work) {
            Ok(h) => hits.extend(h),
            Err(e) => missing.push(format!("AMRFinderPlus failed ({e})")),
        },
        None => missing.push("AMRFinderPlus is not installed on this server".into()),
    }
    match &cfg.vfdb {
        Some(db) => match run_vfdb(tools, db, qry_fasta, work) {
            Ok(h) => hits.extend(h),
            Err(e) => missing.push(format!("the VFDB search failed ({e})")),
        },
        None => missing.push("the VFDB gene set is not installed on this server".into()),
    }
    hits.sort_by(|a, b| (a.contig.as_str(), a.start).cmp(&(b.contig.as_str(), b.start)));
    let status = if missing.len() == 2 {
        ScreenStatus::Unavailable(format!(
            "The resistance and virulence screen did not run: {}.",
            missing.join(" and ")
        ))
    } else if missing.is_empty() {
        ScreenStatus::Done
    } else {
        ScreenStatus::Partial(format!("Only part of the screen ran: {}.", missing[0]))
    };
    Ok((hits, status))
}

/// Run a tool with its own directory first on PATH: AMRFinderPlus calls
/// hmmsearch, blastn and friends by name, installed beside it.
fn command_beside(bin: &Path) -> Command {
    let mut cmd = Command::new(bin);
    if let Some(dir) = bin.parent() {
        let path = std::env::var_os("PATH").unwrap_or_default();
        let mut dirs = vec![dir.to_path_buf()];
        dirs.extend(std::env::split_paths(&path));
        if let Ok(joined) = std::env::join_paths(dirs) {
            cmd.env("PATH", joined);
        }
    }
    cmd
}

/// AMRFinderPlus' `--organism` choices, asked once per process.
fn amr_organisms(bin: &Path, cfg: &ScreenConfig) -> &'static [String] {
    static LIST: OnceLock<Vec<String>> = OnceLock::new();
    LIST.get_or_init(|| {
        let mut cmd = command_beside(bin);
        cmd.arg("-l");
        if let Some(d) = &cfg.amr_db {
            cmd.arg("-d").arg(d);
        }
        let out = cmd.output().ok();
        let text = out
            .map(|o| format!("{}{}", String::from_utf8_lossy(&o.stdout), String::from_utf8_lossy(&o.stderr)))
            .unwrap_or_default();
        parse_organism_list(&text)
    })
}

fn parse_organism_list(text: &str) -> Vec<String> {
    text.lines()
        .find_map(|l| l.split_once("--organism options:").map(|(_, r)| r))
        .map(|r| {
            r.split(',')
                .map(|o| o.trim().to_string())
                .filter(|o| !o.is_empty())
                .collect()
        })
        .unwrap_or_default()
}

/// The `--organism` value for a project organism: "Staphylococcus
/// aureus" -> Staphylococcus_aureus, "Salmonella enterica" -> Salmonella.
/// None for organisms AMRFinderPlus has no species rules for (Listeria,
/// Bacillus): it then searches all genes without point-mutation rules.
pub fn amr_organism(organism: &str, choices: &[String]) -> Option<String> {
    let mut words = organism.split_whitespace();
    let genus = words.next()?;
    let species = format!("{genus}_{}", words.next().unwrap_or(""));
    choices
        .iter()
        .find(|c| c.eq_ignore_ascii_case(&species))
        .or_else(|| choices.iter().find(|c| c.eq_ignore_ascii_case(genus)))
        .cloned()
}

fn run_amrfinder(
    bin: &Path,
    cfg: &ScreenConfig,
    qry_fasta: &Path,
    organism: Option<&str>,
    work: &Path,
) -> Result<Vec<ScreenHit>> {
    let out = work.join("amrfinder.tsv");
    let mut cmd = command_beside(bin);
    cmd.arg("-n")
        .arg(qry_fasta)
        .arg("--plus")
        .args(["--threads", &cfg.threads.max(1).to_string()])
        .arg("-o")
        .arg(&out);
    if let Some(d) = &cfg.amr_db {
        cmd.arg("-d").arg(d);
    }
    if let Some(o) = organism.and_then(|o| amr_organism(o, amr_organisms(bin, cfg))) {
        cmd.args(["-O", &o]);
    }
    let o = cmd
        .output()
        .map_err(|e| EngineError::ToolMissing(format!("amrfinder: {e}")))?;
    if !o.status.success() {
        let msg = String::from_utf8_lossy(&o.stderr);
        let last = msg.lines().rev().find(|l| !l.trim().is_empty()).unwrap_or("");
        return Err(friendly(format!(
            "The resistance screen (AMRFinderPlus) failed. {last}"
        )));
    }
    Ok(parse_amrfinder(&String::from_utf8_lossy(&std::fs::read(&out)?)))
}

/// AMRFinderPlus TSV, by column name: version 4 and the older version 3
/// headers both read.
pub fn parse_amrfinder(text: &str) -> Vec<ScreenHit> {
    let mut lines = text.lines();
    let Some(header) = lines.next() else {
        return Vec::new();
    };
    let cols: Vec<&str> = header.split('\t').collect();
    let idx = |names: &[&str]| {
        cols.iter()
            .position(|c| names.iter().any(|n| c.starts_with(n)))
    };
    let i_gene = idx(&["Element symbol", "Gene symbol"]);
    let i_name = idx(&["Element name", "Sequence name"]);
    let i_type = idx(&["Type", "Element type"]);
    let i_sub = idx(&["Subtype", "Element subtype"]);
    let i_class = idx(&["Class"]);
    let i_subclass = idx(&["Subclass"]);
    let i_contig = idx(&["Contig id"]);
    let i_start = idx(&["Start"]);
    let i_stop = idx(&["Stop"]);
    let i_strand = idx(&["Strand"]);
    let i_method = idx(&["Method"]);
    let i_cov = idx(&["% Coverage of reference"]);
    let i_id = idx(&["% Identity to reference"]);
    let i_acc = idx(&["Closest reference accession", "Accession of closest sequence"]);
    let i_ref = idx(&["Closest reference name", "Name of closest sequence"]);
    let mut hits = Vec::new();
    for line in lines {
        let f: Vec<&str> = line.split('\t').collect();
        let get = |i: Option<usize>| i.and_then(|i| f.get(i)).copied().unwrap_or("").trim();
        let num = |i: Option<usize>| get(i).parse::<f64>().unwrap_or(0.0);
        let gene = get(i_gene);
        if gene.is_empty() {
            continue;
        }
        let class = match (get(i_class), get(i_subclass)) {
            (c, s) if !s.is_empty() && s != "NA" && s != c => format!("{c} / {s}"),
            (c, _) if c != "NA" => c.to_string(),
            _ => String::new(),
        };
        hits.push(ScreenHit {
            source: "AMRFinderPlus".into(),
            gene: gene.to_string(),
            product: get(i_name).to_string(),
            kind: get(i_type).to_string(),
            category: get(i_sub).to_string(),
            class,
            contig: get(i_contig).to_string(),
            start: num(i_start) as u64,
            end: num(i_stop) as u64,
            strand: if get(i_strand) == "-" { -1 } else { 1 },
            identity: num(i_id),
            coverage: num(i_cov),
            reference: format!("{} {}", get(i_acc), get(i_ref)).trim().to_string(),
            method: get(i_method).to_string(),
        });
    }
    hits
}

/// One VFDB core-set record header:
/// `>VFG000068(gb|NP_463735) (actA) actin-assembly inducing protein
/// precursor [ActA (VF0066) - Motility (VFC0204)] [Listeria monocytogenes EGD-e]`
/// -> (id, gene, product, factor, category, organism).
pub fn parse_vfdb_header(h: &str) -> (String, String, String, String, String, String) {
    let h = h.trim_start_matches('>');
    let (id, rest) = h.split_once(' ').unwrap_or((h, ""));
    let rest = rest.trim();
    let (gene, rest) = match rest.strip_prefix('(').and_then(|r| r.split_once(')')) {
        Some((g, r)) => (g.to_string(), r.trim()),
        None => (String::new(), rest),
    };
    // the bracketed groups from the end: [factor - category] [organism]
    let mut groups: Vec<&str> = Vec::new();
    let mut product_end = rest.len();
    let mut s = rest;
    while let Some(close) = s.rfind(']') {
        let Some(open) = s[..close].rfind('[') else { break };
        groups.push(&s[open + 1..close]);
        product_end = open;
        s = &s[..open];
        if groups.len() == 2 {
            break;
        }
    }
    let product = rest[..product_end.min(rest.len())].trim().to_string();
    let organism = groups.first().map(|g| g.to_string()).unwrap_or_default();
    let (factor, category) = groups
        .get(1)
        .map(|g| {
            let (f, c) = g.split_once(" - ").unwrap_or((g, ""));
            let strip = |x: &str| x.split(" (VF").next().unwrap_or(x).trim().to_string();
            (strip(f), strip(c))
        })
        .unwrap_or_default();
    (id.to_string(), gene, product, factor, category, organism)
}

fn run_vfdb(
    tools: &ToolPaths,
    vfdb: &Path,
    qry_fasta: &Path,
    work: &Path,
) -> Result<Vec<ScreenHit>> {
    let db = work.join("screen_db");
    let o = Command::new(&tools.makeblastdb)
        .arg("-in")
        .arg(qry_fasta)
        .args(["-dbtype", "nucl", "-out"])
        .arg(&db)
        .output()
        .map_err(|e| EngineError::ToolMissing(format!("makeblastdb: {e}")))?;
    if !o.status.success() {
        return Err(friendly("The virulence screen could not prepare the genome."));
    }
    let out = work.join("vfdb.tsv");
    let o = Command::new(&tools.blastn)
        .arg("-query")
        .arg(vfdb)
        .arg("-db")
        .arg(&db)
        .args([
            "-evalue",
            "1e-20",
            "-outfmt",
            "6 qseqid sseqid pident length qlen sstart send bitscore",
            "-out",
        ])
        .arg(&out)
        .output()
        .map_err(|e| EngineError::ToolMissing(format!("blastn: {e}")))?;
    if !o.status.success() {
        return Err(friendly(format!(
            "The virulence screen (VFDB) failed. {}",
            String::from_utf8_lossy(&o.stderr).trim()
        )));
    }
    // headers by record id, for the names; read leniently, VFDB carries
    // a few Latin-1 characters in organism names
    let headers = vfdb_headers(&std::fs::read(vfdb)?);
    Ok(select_vfdb_hits(&String::from_utf8_lossy(&std::fs::read(&out)?), &headers))
}

/// Record id -> full header line (without '>'), tolerating bytes that
/// are not UTF-8.
pub fn vfdb_headers(bytes: &[u8]) -> HashMap<String, String> {
    String::from_utf8_lossy(bytes)
        .lines()
        .filter_map(|l| l.strip_prefix('>'))
        .map(|h| {
            let id = h.split_whitespace().next().unwrap_or("").to_string();
            (id, h.trim().to_string())
        })
        .collect()
}

/// VFDB calls from blast rows: identity and coverage of the VFDB gene at
/// least 80 %, and where several VFDB genes hit the same stretch of the
/// genome (homologues from different species), only the best one.
pub fn select_vfdb_hits(text: &str, headers: &HashMap<String, String>) -> Vec<ScreenHit> {
    struct Row {
        q: String,
        s: String,
        pid: f64,
        cov: f64,
        lo: u64,
        hi: u64,
        rev: bool,
        bits: f64,
    }
    let mut rows: Vec<Row> = text
        .lines()
        .filter_map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            if f.len() < 8 {
                return None;
            }
            let (ss, se): (u64, u64) = (f[5].parse().ok()?, f[6].parse().ok()?);
            let len: f64 = f[3].parse().ok()?;
            let qlen: f64 = f[4].parse().ok()?;
            Some(Row {
                q: f[0].to_string(),
                s: f[1].to_string(),
                pid: f[2].parse().ok()?,
                cov: 100.0 * len / qlen.max(1.0),
                lo: ss.min(se),
                hi: ss.max(se),
                rev: se < ss,
                bits: f[7].parse().ok()?,
            })
        })
        .filter(|r| r.pid >= VFDB_MIN_ID && r.cov >= VFDB_MIN_COV)
        .collect();
    rows.sort_by(|a, b| b.bits.partial_cmp(&a.bits).unwrap_or(std::cmp::Ordering::Equal));
    let mut kept: Vec<Row> = Vec::new();
    for r in rows {
        let overlaps = kept.iter().any(|k| {
            k.s == r.s && {
                let ov = k.hi.min(r.hi).saturating_sub(k.lo.max(r.lo));
                ov as f64 > 0.5 * (r.hi - r.lo + 1) as f64
            }
        });
        if !overlaps {
            kept.push(r);
        }
    }
    kept.into_iter()
        .map(|r| {
            let h = headers.get(&r.q).cloned().unwrap_or_else(|| r.q.clone());
            let (id, gene, product, factor, category, organism) = parse_vfdb_header(&h);
            ScreenHit {
                source: "VFDB".into(),
                gene: if gene.is_empty() { id.clone() } else { gene },
                product,
                kind: "VIRULENCE".into(),
                category,
                class: factor,
                contig: r.s,
                start: r.lo,
                end: r.hi,
                strand: if r.rev { -1 } else { 1 },
                identity: r.pid,
                coverage: r.cov.min(100.0),
                reference: format!("{id} [{organism}]"),
                method: "blastn".into(),
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reads_amrfinder_v4_output() {
        let tsv = "Protein id\tContig id\tStart\tStop\tStrand\tElement symbol\tElement name\tScope\tType\tSubtype\tClass\tSubclass\tMethod\tTarget length\tReference sequence length\t% Coverage of reference\t% Identity to reference\tAlignment length\tClosest reference accession\tClosest reference name\tHMM accession\tHMM description\n\
            NA\tLM226_contig1\t2030\t2413\t-\temrC\tmultidrug efflux transporter outer membrane subunit EmrC\tplus\tAMR\tAMR\tQUATERNARY AMMONIUM\tQUATERNARY AMMONIUM\tEXACTX\t128\t128\t100.00\t100.00\t128\tWP_001\temrC\tNA\tNA\n\
            NA\tLM226_contig46\t2466\t2822\t-\tcadC\tCd(II)-sensing repressor CadC\tplus\tSTRESS\tMETAL\tCADMIUM\tCADMIUM\tBLASTX\t119\t120\t99.17\t98.32\t119\tWP_002\tcadC\tNA\tNA\n";
        let h = parse_amrfinder(tsv);
        assert_eq!(h.len(), 2);
        assert_eq!((h[0].gene.as_str(), h[0].kind.as_str(), h[0].class.as_str()), ("emrC", "AMR", "QUATERNARY AMMONIUM"));
        assert_eq!((h[0].start, h[0].end, h[0].strand), (2030, 2413, -1));
        assert_eq!((h[1].category.as_str(), h[1].identity, h[1].coverage), ("METAL", 98.32, 99.17));
        assert!(parse_amrfinder("").is_empty());
    }

    #[test]
    fn picks_amrfinder_organism_only_when_it_has_rules() {
        let list = parse_organism_list("Available --organism options: Campylobacter, Escherichia, Salmonella, Staphylococcus_aureus");
        assert_eq!(list.len(), 4);
        assert_eq!(amr_organism("Staphylococcus aureus", &list).as_deref(), Some("Staphylococcus_aureus"));
        assert_eq!(amr_organism("Salmonella enterica", &list).as_deref(), Some("Salmonella"));
        assert_eq!(amr_organism("Listeria monocytogenes", &list), None);
        assert_eq!(amr_organism("Bacillus cereus", &list), None);
        assert_eq!(amr_organism("", &list), None);
    }

    #[test]
    fn reads_vfdb_headers() {
        let (id, gene, product, factor, cat, org) = parse_vfdb_header(
            ">VFG000068(gb|NP_463735) (actA) actin-assembly inducing protein precursor [ActA (VF0066) - Motility (VFC0204)] [Listeria monocytogenes EGD-e]",
        );
        assert_eq!(id, "VFG000068(gb|NP_463735)");
        assert_eq!(gene, "actA");
        assert_eq!(product, "actin-assembly inducing protein precursor");
        assert_eq!((factor.as_str(), cat.as_str()), ("ActA", "Motility"));
        assert_eq!(org, "Listeria monocytogenes EGD-e");
        let (_, g, p, f, c, _) = parse_vfdb_header("VFG1 (hbp1/svpA) Haemoglobin binding protein 1 [SvpA (VF0263) - Nutritional/Metabolic factor (VFC0272)] [Listeria monocytogenes EGD-e]");
        assert_eq!((g.as_str(), p.as_str(), f.as_str(), c.as_str()), ("hbp1/svpA", "Haemoglobin binding protein 1", "SvpA", "Nutritional/Metabolic factor"));
    }

    #[test]
    fn reads_vfdb_headers_that_are_not_utf8() {
        let mut bytes = b">VFG1(gb|X) (abc) protein [F (VF1) - Cat (VFC1)] [Strain M".to_vec();
        bytes.push(0xfc); // Latin-1 u-umlaut
        bytes.extend_from_slice(b"ller]\nACGT\n>VFG2 (d) q [G (VF2) - C2 (VFC2)] [Y]\nACGT\n");
        let h = vfdb_headers(&bytes);
        assert_eq!(h.len(), 2);
        assert!(h["VFG1(gb|X)"].starts_with("VFG1(gb|X) (abc) protein"));
        assert_eq!(parse_vfdb_header(&h["VFG2"]).1, "d");
    }

    #[test]
    fn keeps_the_best_vfdb_gene_per_genome_stretch() {
        let mut h = HashMap::new();
        h.insert("A".to_string(), "A (hly) listeriolysin O [LLO (VF0064) - Exotoxin (VFC0235)] [Listeria monocytogenes EGD-e]".to_string());
        h.insert("B".to_string(), "B (hly) listeriolysin O [LLO (VF0064) - Exotoxin (VFC0235)] [Listeria ivanovii]".to_string());
        h.insert("C".to_string(), "C (prfA) regulator [PrfA (VF0062) - Regulation (VFC0301)] [Listeria monocytogenes EGD-e]".to_string());
        let rows = "A\tc1\t99.9\t1590\t1590\t100\t1689\t2900\n\
            B\tc1\t85.0\t1580\t1590\t105\t1684\t1800\n\
            C\tc1\t70.0\t714\t714\t5000\t5713\t500\n\
            C\tc2\t99.0\t300\t714\t1\t300\t500\n";
        let hits = select_vfdb_hits(rows, &h);
        // B overlaps A (homologue), C fails identity on c1 and coverage on c2
        assert_eq!(hits.len(), 1);
        assert_eq!((hits[0].gene.as_str(), hits[0].class.as_str(), hits[0].category.as_str()), ("hly", "LLO", "Exotoxin"));
        assert_eq!(hits[0].reference, "A [Listeria monocytogenes EGD-e]");
    }
}
