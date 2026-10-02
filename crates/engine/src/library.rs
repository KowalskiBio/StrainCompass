//! The local reference library of a genus: its complete RefSeq genomes,
//! near-duplicates removed, with every gene sorted into variant groups and
//! indexed for BLAST. Built offline by `tools/library/build_library.py`;
//! the app only reads it.
//!
//! It answers what a gene name alone cannot: which of the genus' several
//! genes of that name exist, where each sits (which plasmid, the
//! chromosome) and how common it is. "cadA" in Listeria is the chromosomal
//! CadA, the cadA of Tn5422 and a pLI100-type one; a panel holds them all
//! as variants and the search reports which one the strains carry.
//!
//! Installed as `<root>/<genus, lower case>/current`, a link to one
//! version folder, so a new build sits beside the old one until the link
//! is switched. Optional, like the resistance screen: no library means
//! the app falls back to the curated catalogs and NCBI.

use crate::{friendly, Result};
use rusqlite::{Connection, OpenFlags, OptionalExtension};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::process::Command;

/// At most this many variants are offered for one name; the most common
/// first. A name rarely covers more than three or four real variants.
pub const MAX_VARIANTS: usize = 8;

/// An unnamed group resembling a named one is offered as a variant of the
/// name from this protein identity on. Any variant matching calls the
/// gene present, so the bar is high: the pLI100-type cadA (76 % to the
/// cadA of Tn5422) passes, the internalins resembling inlA (50-68 %) do
/// not. The library keeps the weaker resemblances (from 50 %) as evidence.
pub const MIN_HOMOLOG_IDENTITY: f64 = 70.0;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FileEntry {
    pub sha256: String,
    pub bytes: u64,
}

/// `manifest.json`: what was built, from what, how.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Manifest {
    pub format: u32,
    pub genus: String,
    pub version: String,
    pub built: String,
    #[serde(default)]
    pub amrfinder_db: String,
    #[serde(default)]
    pub counts: BTreeMap<String, u64>,
    #[serde(default)]
    pub files: BTreeMap<String, FileEntry>,
}

pub struct Library {
    pub dir: PathBuf,
    pub manifest: Manifest,
}

/// One variant of a gene name: a variant group of the library.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct LibraryVariant {
    pub group: i64,
    /// The group's own gene name, or the requested one when it has none.
    pub label: String,
    pub product: String,
    /// How the name led here: "gene" (annotated so), "locus" (locus tag),
    /// "catalog" (a curated entry or seed record matched its protein),
    /// "homolog" (unnamed, resembling a group of that name).
    pub matched_by: String,
    /// Protein identity behind the match; 100 for annotated names.
    pub identity: f64,
    /// What made the match, in words ("GenBank L28104.1 cadA",
    /// "76 % identical to the cadA of group 1047").
    pub evidence: String,
    /// The group's representative gene, DNA, in its own orientation.
    pub seq: Vec<u8>,
    /// Where the representative sits: "NZ_CP013725.1:47228-49363".
    pub locus: String,
    /// Its replicon in words: "plasmid pLmN1546 (Listeria monocytogenes 2015TE24968)".
    pub place: String,
    /// Genomes carrying the gene, counting the near-duplicates each kept
    /// replicon stands for.
    pub n_genomes: u64,
    pub n_plasmid: u64,
    pub n_chromosome: u64,
}

/// A replicon matching a sequence.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct LibraryHit {
    pub accession: String,
    /// "chromosome" or "plasmid".
    pub kind: String,
    /// The plasmid name, empty for a chromosome or an unnamed plasmid.
    pub name: String,
    pub organism: String,
    pub length: u64,
    pub identity: f64,
    /// Share of the query the best match covers, percent.
    pub coverage: f64,
    pub start: u64,
    pub end: u64,
    /// The near-identical replicons this one stands for, itself included.
    pub represents: u64,
}

impl Library {
    /// The installed library of `genus` under `root`, None when there is
    /// none; an error when one is there but broken.
    pub fn for_genus(root: &Path, genus: &str) -> Option<Result<Library>> {
        let genus = genus.trim().to_ascii_lowercase();
        if genus.is_empty() {
            return None;
        }
        let dir = root.join(&genus).join("current");
        dir.join("manifest.json")
            .is_file()
            .then(|| Library::open(&dir))
    }

    /// Open a library folder. Every file the manifest lists must be there
    /// at its listed size, so a half-copied folder is refused rather than
    /// giving wrong answers.
    pub fn open(dir: &Path) -> Result<Library> {
        let text = std::fs::read_to_string(dir.join("manifest.json")).map_err(|e| {
            friendly(format!(
                "The reference library at {} has no readable manifest. ({e})",
                dir.display()
            ))
        })?;
        let manifest: Manifest = serde_json::from_str(&text).map_err(|e| {
            friendly(format!(
                "The reference library manifest is not valid. ({e})"
            ))
        })?;
        if manifest.format != 1 {
            return Err(friendly(format!(
                "The reference library at {} has format {}, this app reads format 1.",
                dir.display(),
                manifest.format
            )));
        }
        for (name, f) in &manifest.files {
            let len = std::fs::metadata(dir.join(name)).map(|m| m.len()).ok();
            if len != Some(f.bytes) {
                return Err(friendly(format!(
                    "The reference library at {} is incomplete: {name} is missing or has the wrong size. Copy the folder again.",
                    dir.display()
                )));
            }
        }
        Ok(Library {
            dir: dir.to_path_buf(),
            manifest,
        })
    }

    /// "Listeria library 2026-10-02"
    pub fn title(&self) -> String {
        format!("{} library {}", self.manifest.genus, self.manifest.version)
    }

    fn conn(&self) -> Result<Connection> {
        Connection::open_with_flags(
            self.dir.join("library.sqlite"),
            OpenFlags::SQLITE_OPEN_READ_ONLY | OpenFlags::SQLITE_OPEN_NO_MUTEX,
        )
        .map_err(sql)
    }

    /// The variants of a gene name (symbol, locus tag or old locus tag such
    /// as lmo0444): named matches first, then resembling unnamed groups,
    /// each the most common first.
    pub fn variants(&self, name: &str) -> Result<Vec<LibraryVariant>> {
        let name = name.trim();
        if name.is_empty() {
            return Ok(Vec::new());
        }
        let conn = self.conn()?;
        // one row per group: its strongest kind of match
        let mut stmt = conn
            .prepare(
                "SELECT n.group_id, n.kind, n.identity, n.via_group, n.name
                 FROM names n WHERE n.name = ?1 COLLATE NOCASE",
            )
            .map_err(sql)?;
        let rows = stmt
            .query_map([name], |r| {
                Ok((
                    r.get::<_, i64>(0)?,
                    r.get::<_, String>(1)?,
                    r.get::<_, Option<f64>>(2)?.unwrap_or(0.0),
                    r.get::<_, Option<i64>>(3)?,
                    r.get::<_, String>(4)?,
                ))
            })
            .map_err(sql)?;
        // the library's spelling of the name ("cadA" for "CADA")
        let mut spelled = name.to_string();
        let mut best: BTreeMap<i64, (u8, String, f64, Option<i64>)> = BTreeMap::new();
        for row in rows {
            let (group, kind, identity, via, stored) = row.map_err(sql)?;
            if kind == "homolog" && identity < MIN_HOMOLOG_IDENTITY {
                continue;
            }
            if stored == name || spelled == name {
                spelled = stored;
            }
            let rank = kind_rank(&kind);
            let better = best
                .get(&group)
                .is_none_or(|(r, _, i, _)| rank < *r || (rank == *r && identity > *i));
            if better {
                best.insert(group, (rank, kind, identity, via));
            }
        }
        let mut out = Vec::new();
        for (group, (_, kind, identity, via)) in best {
            out.push(self.variant(&conn, group, &spelled, &kind, identity, via)?);
        }
        out.sort_by(|a, b| {
            kind_rank(&a.matched_by)
                .cmp(&kind_rank(&b.matched_by))
                .then(b.n_genomes.cmp(&a.n_genomes))
                .then(a.group.cmp(&b.group))
        });
        out.truncate(MAX_VARIANTS);
        Ok(out)
    }

    fn variant(
        &self,
        conn: &Connection,
        group: i64,
        name: &str,
        kind: &str,
        identity: f64,
        via: Option<i64>,
    ) -> Result<LibraryVariant> {
        let (label, product, nt, rep_gene, n_genomes, n_plasmid, n_chromosome) = conn
            .query_row(
                "SELECT label, product, nt, rep_gene, n_genomes, n_plasmid, n_chromosome
                 FROM groups WHERE id = ?1",
                [group],
                |r| {
                    Ok((
                        r.get::<_, String>(0)?,
                        r.get::<_, String>(1)?,
                        r.get::<_, String>(2)?,
                        r.get::<_, i64>(3)?,
                        r.get::<_, i64>(4)? as u64,
                        r.get::<_, i64>(5)? as u64,
                        r.get::<_, i64>(6)? as u64,
                    ))
                },
            )
            .map_err(sql)?;
        let (locus, place) = conn
            .query_row(
                "SELECT g.replicon, g.start, g.end, r.kind, r.name, a.organism, a.strain
                 FROM genes g JOIN replicons r ON r.accession = g.replicon
                 LEFT JOIN assemblies a ON a.accession = r.assembly
                 WHERE g.id = ?1",
                [rep_gene],
                |r| {
                    let acc: String = r.get(0)?;
                    let (s, e): (i64, i64) = (r.get(1)?, r.get(2)?);
                    let kind: String = r.get(3)?;
                    let pname: String = r.get(4)?;
                    let org: Option<String> = r.get(5)?;
                    let strain: Option<String> = r.get(6)?;
                    Ok((
                        format!("{acc}:{s}-{e}"),
                        replicon_words(&kind, &pname, &with_strain(org, strain)),
                    ))
                },
            )
            .map_err(sql)?;
        let evidence = match kind {
            "catalog" => {
                let hit: Option<(String, String, f64)> = conn
                    .query_row(
                        "SELECT source, symbol, identity FROM catalog_hits
                         WHERE group_id = ?1 ORDER BY identity DESC LIMIT 1",
                        [group],
                        |r| Ok((r.get(0)?, r.get(1)?, r.get(2)?)),
                    )
                    .optional()
                    .map_err(sql)?;
                match hit {
                    Some((src, sym, id)) => {
                        format!("matches {src} {sym} ({id:.0} % protein identity)")
                    }
                    None => "matches a curated entry".into(),
                }
            }
            "homolog" => {
                let other = via
                    .and_then(|v| {
                        conn.query_row("SELECT label FROM groups WHERE id = ?1", [v], |r| {
                            r.get::<_, String>(0)
                        })
                        .ok()
                    })
                    .filter(|l| !l.is_empty())
                    .unwrap_or_else(|| name.to_string());
                format!(
                    "not named in its annotation; {identity:.0} % protein identity to the {other} of group {}",
                    via.unwrap_or(0)
                )
            }
            "locus" => format!("locus tag {name}"),
            _ => "annotated with this name".into(),
        };
        Ok(LibraryVariant {
            group,
            label: if label.is_empty() {
                name.to_string()
            } else {
                label
            },
            product,
            matched_by: kind.to_string(),
            identity,
            evidence,
            seq: nt.into_bytes(),
            locus,
            place,
            n_genomes,
            n_plasmid,
            n_chromosome,
        })
    }

    /// The replicons carrying `seq`, best match each, strongest first.
    /// Runs blastn against the library's genomes in `work`.
    pub fn search(&self, blastn: &Path, seq: &[u8], work: &Path) -> Result<Vec<LibraryHit>> {
        std::fs::create_dir_all(work)?;
        let query = work.join("library_query.fasta");
        let mut fasta = b">query\n".to_vec();
        fasta.extend_from_slice(seq);
        fasta.push(b'\n');
        std::fs::write(&query, fasta)?;
        let out = Command::new(blastn)
            .arg("-query")
            .arg(&query)
            .arg("-db")
            .arg(self.dir.join("blast").join("genomes"))
            .args([
                "-outfmt",
                "6 sseqid pident length qlen sstart send bitscore",
                "-evalue",
                "1e-20",
                "-max_target_seqs",
                "5000",
                "-num_threads",
                "2",
            ])
            .output()
            .map_err(|e| crate::EngineError::ToolMissing(format!("blastn: {e}")))?;
        if !out.status.success() {
            return Err(friendly(format!(
                "The reference library search failed. {}",
                String::from_utf8_lossy(&out.stderr).trim()
            )));
        }
        // best-scoring match per replicon
        let mut best: BTreeMap<String, (f64, f64, f64, u64, u64)> = BTreeMap::new();
        for line in String::from_utf8_lossy(&out.stdout).lines() {
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 7 {
                continue;
            }
            let acc = f[0]
                .trim_start_matches("ref|")
                .trim_end_matches('|')
                .to_string();
            let num = |i: usize| f[i].parse::<f64>().unwrap_or(0.0);
            let (pid, len, qlen, bits) = (num(1), num(2), num(3), num(6));
            let (s, e) = (num(4) as u64, num(5) as u64);
            let cov = if qlen > 0.0 {
                (100.0 * len / qlen).min(100.0)
            } else {
                0.0
            };
            if best.get(&acc).is_none_or(|b| bits > b.0) {
                best.insert(acc, (bits, pid, cov, s.min(e), s.max(e)));
            }
        }
        let conn = self.conn()?;
        let mut hits = Vec::new();
        for (acc, (_, identity, coverage, start, end)) in best {
            let row = conn
                .query_row(
                    "SELECT r.kind, r.name, r.length, r.represents, a.organism, a.strain
                     FROM replicons r LEFT JOIN assemblies a ON a.accession = r.assembly
                     WHERE r.accession = ?1",
                    [&acc],
                    |r| {
                        Ok((
                            r.get::<_, String>(0)?,
                            r.get::<_, String>(1)?,
                            r.get::<_, i64>(2)? as u64,
                            r.get::<_, i64>(3)? as u64,
                            with_strain(r.get(4)?, r.get(5)?),
                        ))
                    },
                )
                .optional()
                .map_err(sql)?;
            let Some((kind, name, length, represents, organism)) = row else {
                continue;
            };
            hits.push(LibraryHit {
                accession: acc,
                kind,
                name,
                organism,
                length,
                identity,
                coverage,
                start,
                end,
                represents: represents.max(1),
            });
        }
        hits.sort_by(|a, b| {
            (b.coverage * b.identity)
                .total_cmp(&(a.coverage * a.identity))
                .then(a.length.cmp(&b.length))
        });
        Ok(hits)
    }

    /// The variant groups each sequence belongs to: groups whose
    /// representative gene it matches at DNA level over most of both
    /// lengths. A reference gene's own group is not another variant of it.
    pub fn groups_of(
        &self,
        blastn: &Path,
        seqs: &[(String, Vec<u8>)],
        work: &Path,
    ) -> Result<BTreeMap<String, Vec<i64>>> {
        let mut out: BTreeMap<String, Vec<i64>> = BTreeMap::new();
        if seqs.is_empty() {
            return Ok(out);
        }
        std::fs::create_dir_all(work)?;
        let query = work.join("library_groups_query.fasta");
        let mut fasta = Vec::new();
        for (i, (_, seq)) in seqs.iter().enumerate() {
            fasta.extend_from_slice(format!(">q{i}\n").as_bytes());
            fasta.extend_from_slice(seq);
            fasta.push(b'\n');
        }
        std::fs::write(&query, fasta)?;
        let res = Command::new(blastn)
            .arg("-query")
            .arg(&query)
            .arg("-db")
            .arg(self.dir.join("blast").join("groups_nt"))
            .args([
                "-task",
                "dc-megablast",
                "-outfmt",
                "6 qseqid sseqid pident length qlen slen",
                "-evalue",
                "1e-20",
                "-max_target_seqs",
                "50",
            ])
            .output()
            .map_err(|e| crate::EngineError::ToolMissing(format!("blastn: {e}")))?;
        if !res.status.success() {
            return Err(friendly(format!(
                "The reference library search failed. {}",
                String::from_utf8_lossy(&res.stderr).trim()
            )));
        }
        for line in String::from_utf8_lossy(&res.stdout).lines() {
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 6 {
                continue;
            }
            let num = |i: usize| f[i].parse::<f64>().unwrap_or(0.0);
            let cov = 100.0 * num(3) / num(4).max(num(5)).max(1.0);
            if num(2) < GROUP_DNA_IDENTITY || cov < GROUP_DNA_COVERAGE {
                continue;
            }
            let (Some(qi), Some(group)) = (
                f[0].strip_prefix('q').and_then(|n| n.parse::<usize>().ok()),
                group_id(f[1]),
            ) else {
                continue;
            };
            let Some((name, _)) = seqs.get(qi) else {
                continue;
            };
            let held = out.entry(name.clone()).or_default();
            if !held.contains(&group) {
                held.push(group);
            }
        }
        Ok(out)
    }

    /// The protein BLAST database of the variant groups (one protein per
    /// group, titled "g<id> <name or -> <product>").
    pub fn groups_protein_db(&self) -> PathBuf {
        self.dir.join("blast").join("groups")
    }

    /// A variant group in brief: (representative protein accession, gene
    /// name, product, genomes carrying it); None for an unknown id.
    pub fn group_brief(&self, group: i64) -> Result<Option<(String, String, String, u64)>> {
        self.conn()?
            .query_row(
                "SELECT rep_protein, label, product, n_genomes FROM groups WHERE id = ?1",
                [group],
                |r| Ok((r.get(0)?, r.get(1)?, r.get(2)?, r.get::<_, i64>(3)? as u64)),
            )
            .optional()
            .map_err(sql)
    }

    /// A genome stored in the library, by its RefSeq (GCF_) or GenBank
    /// (GCA_) assembly accession: every assembly of the build, kept or
    /// merged, with the genome and annotation files NCBI serves for it.
    pub fn reference(&self, accession: &str) -> Result<Option<StoredReference>> {
        let acc = accession.trim();
        if acc.is_empty() {
            return Ok(None);
        }
        let conn = self.conn()?;
        let row = conn
            .query_row(
                "SELECT accession, assembly_name, organism, strain FROM assemblies
                 WHERE accession = ?1 OR paired_accession = ?1",
                [acc],
                |r| {
                    Ok((
                        r.get::<_, String>(0)?,
                        r.get::<_, Option<String>>(1)?.unwrap_or_default(),
                        r.get::<_, Option<String>>(2)?,
                        r.get::<_, Option<String>>(3)?,
                    ))
                },
            )
            .optional()
            .map_err(sql)?;
        let Some((accession, assembly_name, organism, strain)) = row else {
            return Ok(None);
        };
        let fna_gz = self.reference_file(&accession, "fna");
        let gff_gz = self.reference_file(&accession, "gff");
        if !fna_gz.is_file() || !gff_gz.is_file() {
            return Ok(None);
        }
        Ok(Some(StoredReference {
            accession,
            assembly_name,
            organism: with_strain(organism, strain),
            fna_gz,
            gff_gz,
        }))
    }

    fn reference_file(&self, assembly: &str, ext: &str) -> PathBuf {
        self.dir
            .join("references")
            .join(format!("{assembly}.{ext}.gz"))
    }

    /// The library's copy of a sequence record: a RefSeq accession as is,
    /// a GenBank one through its RefSeq twin (CP013725.1 is
    /// NZ_CP013725.1). Every replicon of every assembly, kept or merged.
    pub fn replicon(&self, accession: &str) -> Result<Option<Replicon>> {
        let acc = accession.trim();
        if acc.is_empty() {
            return Ok(None);
        }
        let conn = self.conn()?;
        for candidate in [acc.to_string(), format!("NZ_{acc}")] {
            let row = conn
                .query_row(
                    "SELECT r.accession, r.assembly, r.kind, r.name, r.length, a.organism, a.strain
                     FROM replicons r LEFT JOIN assemblies a ON a.accession = r.assembly
                     WHERE r.accession = ?1",
                    [&candidate],
                    |r| {
                        Ok(Replicon {
                            accession: r.get(0)?,
                            assembly: r.get(1)?,
                            kind: r.get(2)?,
                            name: r.get(3)?,
                            length: r.get::<_, i64>(4)? as u64,
                            organism: with_strain(r.get(5)?, r.get(6)?),
                        })
                    },
                )
                .optional()
                .map_err(sql)?;
            if row.is_some() {
                return Ok(row);
            }
        }
        Ok(None)
    }

    /// A replicon's sequence, from its assembly's stored genome.
    pub fn replicon_seq(&self, rep: &Replicon) -> Result<Option<Vec<u8>>> {
        let path = self.reference_file(&rep.assembly, "fna");
        if !path.is_file() {
            return Ok(None);
        }
        let text = gunzip_text(&path)?;
        let recs = crate::fasta::parse_fasta_str(&text)?;
        Ok(recs
            .into_iter()
            .find(|r| r.id == rep.accession)
            .map(|r| r.seq))
    }

    /// The annotated genes of a replicon, from its assembly's stored
    /// annotation.
    pub fn replicon_genes(&self, rep: &Replicon) -> Result<Vec<crate::gff::Gene>> {
        let path = self.reference_file(&rep.assembly, "gff");
        if !path.is_file() {
            return Ok(Vec::new());
        }
        let genes = crate::gff::parse_gff_str(&gunzip_text(&path)?)?;
        Ok(genes
            .into_iter()
            .filter(|g| g.seqid == rep.accession)
            .collect())
    }
}

/// A genome of the library, ready to serve as a project reference.
#[derive(Debug, Clone)]
pub struct StoredReference {
    /// The RefSeq assembly accession (GCF_...).
    pub accession: String,
    pub assembly_name: String,
    pub organism: String,
    pub fna_gz: PathBuf,
    pub gff_gz: PathBuf,
}

/// One chromosome or plasmid of the library.
#[derive(Debug, Clone, PartialEq)]
pub struct Replicon {
    /// Its RefSeq accession, even when asked for by the GenBank one.
    pub accession: String,
    pub assembly: String,
    pub kind: String,
    pub name: String,
    pub length: u64,
    pub organism: String,
}

impl Replicon {
    /// "Listeria monocytogenes X plasmid pLM80"
    pub fn title(&self) -> String {
        let what = match (self.kind.as_str(), self.name.is_empty()) {
            ("plasmid", false) => format!("plasmid {}", self.name),
            ("plasmid", true) => "plasmid".to_string(),
            _ => "chromosome".to_string(),
        };
        format!("{} {what}", self.organism).trim().to_string()
    }

    /// `seq[lo..=hi]` (1-based) of this replicon, reverse-complemented for
    /// the minus strand; None when the range does not fit.
    pub fn slice(seq: &[u8], lo: u64, hi: u64, minus: bool) -> Option<Vec<u8>> {
        if lo == 0 || hi < lo || hi as usize > seq.len() {
            return None;
        }
        let s = seq[(lo - 1) as usize..hi as usize].to_ascii_uppercase();
        Some(if minus { crate::fasta::revcomp(&s) } else { s })
    }
}

fn gunzip_text(path: &Path) -> Result<String> {
    use std::io::Read;
    let mut text = String::new();
    flate2::read::GzDecoder::new(std::fs::File::open(path)?).read_to_string(&mut text)?;
    Ok(text)
}

/// A sequence belongs to a variant group when it matches the group's
/// representative gene this closely (groups are 90 % protein identity
/// clusters; their DNA is rarely less alike than this).
const GROUP_DNA_IDENTITY: f64 = 85.0;
const GROUP_DNA_COVERAGE: f64 = 80.0;

/// "g123" (or "lcl|g123") -> 123. BLAST reports parsed local ids in
/// upper case ("G123") for real-size databases.
pub fn group_id(sseqid: &str) -> Option<i64> {
    sseqid
        .rsplit('|')
        .find(|p| !p.is_empty())?
        .strip_prefix(['g', 'G'])?
        .parse()
        .ok()
}

/// Named matches before resembling ones.
fn kind_rank(kind: &str) -> u8 {
    match kind {
        "gene" | "catalog" | "locus" => 0,
        _ => 1,
    }
}

/// "Listeria monocytogenes" + "J1-220" -> "Listeria monocytogenes J1-220";
/// NCBI often already puts the strain in the organism name.
fn with_strain(organism: Option<String>, strain: Option<String>) -> String {
    let organism = organism.unwrap_or_default();
    match strain.filter(|s| !s.is_empty() && !organism.contains(s.as_str())) {
        Some(s) => format!("{organism} {s}").trim().to_string(),
        None => organism,
    }
}

fn replicon_words(kind: &str, name: &str, organism: &str) -> String {
    let what = match (kind, name.is_empty()) {
        ("plasmid", false) => format!("plasmid {name}"),
        ("plasmid", true) => "an unnamed plasmid".to_string(),
        _ => "chromosome".to_string(),
    };
    if organism.is_empty() {
        what
    } else {
        format!("{what} ({organism})")
    }
}

fn sql(e: rusqlite::Error) -> crate::EngineError {
    friendly(format!("The reference library could not be read. ({e})"))
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    /// The tables of `build_library.py`, filled with a tiny cadA world:
    /// the chromosomal CadA (annotated), Tn5422 cadA (named by a seed
    /// record) and a pLI100-type one (unnamed, 76 % to Tn5422's).
    pub fn fixture(dir: &Path) -> Library {
        std::fs::create_dir_all(dir.join("blast")).unwrap();
        let conn = Connection::open(dir.join("library.sqlite")).unwrap();
        conn.execute_batch(
            "CREATE TABLE meta(key TEXT PRIMARY KEY, value TEXT);
             CREATE TABLE assemblies(accession TEXT PRIMARY KEY, organism TEXT, strain TEXT, kept INTEGER, represented_by TEXT, assembly_name TEXT, paired_accession TEXT);
             CREATE TABLE replicons(accession TEXT PRIMARY KEY, assembly TEXT, kind TEXT, name TEXT, length INTEGER, circular INTEGER, kept INTEGER, represented_by TEXT, represents INTEGER);
             CREATE TABLE groups(id INTEGER PRIMARY KEY, label TEXT, product TEXT, rep_protein TEXT, rep_gene INTEGER, protein TEXT, nt TEXT, n_genes INTEGER, n_genomes INTEGER, n_plasmid INTEGER, n_chromosome INTEGER);
             CREATE TABLE genes(id INTEGER PRIMARY KEY, replicon TEXT, start INTEGER, end INTEGER, strand TEXT, name TEXT, locus_tag TEXT, old_locus_tag TEXT, product TEXT, protein_id TEXT, pseudo INTEGER, group_id INTEGER, nt TEXT);
             CREATE TABLE names(name TEXT COLLATE NOCASE, group_id INTEGER, kind TEXT, n INTEGER, identity REAL, via_group INTEGER);
             CREATE TABLE catalog_hits(group_id INTEGER, source TEXT, symbol TEXT, product TEXT, identity REAL, coverage REAL);
             INSERT INTO assemblies VALUES ('GCF_1','Listeria monocytogenes Scott A','Scott A',1,'GCF_1','ASM1v1','GCA_1'),
                                           ('GCF_2','Listeria monocytogenes','X',1,'GCF_2','ASM2v1','GCA_2');
             INSERT INTO replicons VALUES ('NZ_CP1.1','GCF_1','chromosome','',3000000,1,1,'NZ_CP1.1',3),
                                          ('NZ_CP2.1','GCF_2','plasmid','pLmN1546',60000,1,1,'NZ_CP2.1',4),
                                          ('NC_3.1','GCF_2','plasmid','',80000,1,1,'NC_3.1',6);
             INSERT INTO genes VALUES (1,'NZ_CP1.1',100,2199,'+','cadA','T_1','','cadmium-translocating P-type ATPase CadA','WP_1',0,10,NULL),
                                      (2,'NZ_CP2.1',47228,49363,'-','','T_2','','heavy metal translocating P-type ATPase','WP_2',0,11,NULL),
                                      (3,'NC_3.1',61300,63417,'+','','T_3','lmo9999','heavy metal translocating P-type ATPase','WP_3',0,12,NULL);
             INSERT INTO groups VALUES (10,'cadA','cadmium-translocating P-type ATPase CadA','WP_1',1,'M','ATGAAA',1,3,0,3),
                                       (11,'','heavy metal translocating P-type ATPase','WP_2',2,'M','ATGCCC',1,4,4,0),
                                       (12,'','heavy metal translocating P-type ATPase','WP_3',3,'M','ATGGGG',1,6,6,0);
             INSERT INTO names VALUES ('cadA',10,'gene',1,100,NULL),('cadA',10,'catalog',1,100,NULL),
                                      ('cadA',11,'catalog',1,100,NULL),
                                      ('cadA',12,'homolog',1,75.7,11),('lmo9999',12,'locus',1,100,NULL),
                                      ('inlA',10,'homolog',1,67.9,99);
             INSERT INTO catalog_hits VALUES (10,'AMRFinderPlus','cadA_Lm','CadA',100,99.9),
                                             (11,'GenBank L28104.1','cadA','ATPase',100,100);",
        )
        .unwrap();
        drop(conn);
        // GCF_2's stored genome: the pLmN1546 plasmid, cadA on its minus
        // strand at 11-16
        std::fs::create_dir_all(dir.join("references")).unwrap();
        let gz = |name: &str, text: &str| {
            use std::io::Write;
            let f = std::fs::File::create(dir.join("references").join(name)).unwrap();
            let mut e = flate2::write::GzEncoder::new(f, flate2::Compression::fast());
            e.write_all(text.as_bytes()).unwrap();
            e.finish().unwrap();
        };
        gz("GCF_2.fna.gz", ">NZ_CP2.1 Listeria monocytogenes X plasmid pLmN1546\nAAAAAAAAAAGGGCATAAAAAAAA\n>NC_3.1 other\nCCCC\n");
        gz(
            "GCF_2.gff.gz",
            "##gff-version 3\n\
             NZ_CP2.1\tRefSeq\tgene\t11\t16\t.\t-\t.\tID=gene-T_2;gene=cadA;locus_tag=T_2;gene_biotype=protein_coding\n\
             NZ_CP2.1\tRefSeq\tCDS\t11\t16\t.\t-\t0\tID=cds-WP_2;Parent=gene-T_2;locus_tag=T_2;product=ATPase;protein_id=WP_2\n\
             NC_3.1\tRefSeq\tgene\t1\t4\t.\t+\t.\tID=gene-T_3;locus_tag=T_3;gene_biotype=protein_coding\n",
        );
        let bytes = std::fs::metadata(dir.join("library.sqlite")).unwrap().len();
        let manifest = serde_json::json!({
            "format": 1, "genus": "Listeria", "version": "test", "built": "2026-10-02",
            "counts": {"groups": 3},
            "files": {"library.sqlite": {"sha256": "-", "bytes": bytes}},
        });
        std::fs::write(dir.join("manifest.json"), manifest.to_string()).unwrap();
        Library::open(dir).unwrap()
    }

    #[test]
    fn lists_every_variant_of_a_name_named_first() {
        let tmp = tempdir("variants");
        let lib = fixture(&tmp);
        let v = lib.variants("CADA").unwrap();
        let groups: Vec<i64> = v.iter().map(|x| x.group).collect();
        // named (annotation or seed record) first, most common first; the
        // resembling unnamed one last, though it is the most common
        assert_eq!(groups, [11, 10, 12]);
        assert_eq!(v[0].matched_by, "catalog");
        assert!(
            v[0].evidence.contains("GenBank L28104.1"),
            "{}",
            v[0].evidence
        );
        assert_eq!(v[0].label, "cadA");
        assert_eq!(v[0].place, "plasmid pLmN1546 (Listeria monocytogenes X)");
        assert_eq!(v[0].locus, "NZ_CP2.1:47228-49363");
        assert_eq!(v[1].matched_by, "gene");
        assert_eq!(v[1].place, "chromosome (Listeria monocytogenes Scott A)");
        assert_eq!(v[2].matched_by, "homolog");
        assert_eq!(v[2].identity, 75.7);
        assert!(
            v[2].evidence.contains("76 % protein identity"),
            "{}",
            v[2].evidence
        );
        // a weak resemblance is no variant
        assert!(lib.variants("inlA").unwrap().is_empty());
        assert_eq!(v[2].place, "an unnamed plasmid (Listeria monocytogenes X)");
        assert_eq!(v[2].seq, b"ATGGGG");
    }

    #[test]
    fn finds_old_locus_tags_and_nothing_for_unknown_names() {
        let tmp = tempdir("locus");
        let lib = fixture(&tmp);
        let v = lib.variants("lmo9999").unwrap();
        assert_eq!(v.len(), 1);
        assert_eq!((v[0].group, v[0].matched_by.as_str()), (12, "locus"));
        assert!(lib.variants("nope").unwrap().is_empty());
        assert!(lib.variants("  ").unwrap().is_empty());
    }

    #[test]
    fn refuses_a_half_copied_library() {
        let tmp = tempdir("broken");
        fixture(&tmp);
        std::fs::write(tmp.join("library.sqlite"), b"short").unwrap();
        let err = Library::open(&tmp).err().unwrap().to_string();
        assert!(err.contains("incomplete"), "{err}");
    }

    #[test]
    fn is_found_per_genus_through_the_current_link() {
        let root = tempdir("root");
        assert!(Library::for_genus(&root, "Listeria").is_none());
        let v = root.join("listeria").join("test");
        fixture(&v);
        #[cfg(unix)]
        std::os::unix::fs::symlink("test", root.join("listeria").join("current")).unwrap();
        let lib = Library::for_genus(&root, "Listeria").unwrap().unwrap();
        assert_eq!(lib.title(), "Listeria library test");
        assert!(Library::for_genus(&root, "Bacillus").is_none());
        assert!(Library::for_genus(&root, "").is_none());
    }

    #[test]
    fn serves_stored_genomes_and_records() {
        let tmp = tempdir("stored");
        let lib = fixture(&tmp);
        // by the GenBank twin of the assembly accession
        let r = lib.reference("GCA_2").unwrap().unwrap();
        assert_eq!(
            (r.accession.as_str(), r.assembly_name.as_str()),
            ("GCF_2", "ASM2v1")
        );
        assert!(r.fna_gz.is_file() && r.gff_gz.is_file());
        // a stored assembly without files is not served
        assert!(lib.reference("GCF_1").unwrap().is_none());
        assert!(lib.reference("GCF_9").unwrap().is_none());

        // a GenBank record accession through its RefSeq twin
        let rep = lib.replicon("CP2.1").unwrap().unwrap();
        assert_eq!(rep.accession, "NZ_CP2.1");
        assert_eq!(rep.title(), "Listeria monocytogenes X plasmid pLmN1546");
        let seq = lib.replicon_seq(&rep).unwrap().unwrap();
        assert_eq!(seq.len(), 24);
        let genes = lib.replicon_genes(&rep).unwrap();
        assert_eq!(genes.len(), 1);
        assert_eq!(
            (genes[0].symbol.as_str(), genes[0].product.as_str()),
            ("cadA", "ATPase")
        );
        let g = &genes[0];
        assert_eq!(
            Replicon::slice(&seq, g.start, g.end, g.strand < 0).unwrap(),
            b"ATGCCC"
        );
        assert_eq!(Replicon::slice(&seq, 20, 30, false), None);
        assert!(lib.replicon("XX1.1").unwrap().is_none());
    }

    #[test]
    fn reads_group_ids_from_blast_subjects() {
        assert_eq!(group_id("g123"), Some(123));
        assert_eq!(group_id("lcl|g7"), Some(7));
        assert_eq!(group_id("G12299"), Some(12299));
        assert_eq!(group_id("NZ_CP1.1"), None);
    }

    /// With BLAST+ installed: a sequence finds its own group, and a
    /// replicon search finds the plasmid carrying it.
    #[test]
    fn searches_the_blast_databases_when_blast_is_installed() {
        let (Some(makeblastdb), Some(blastn)) = (
            crate::tools::find_optional("makeblastdb"),
            crate::tools::find_optional("blastn"),
        ) else {
            eprintln!("BLAST+ not installed; skipped");
            return;
        };
        let tmp = tempdir("blast");
        let lib = fixture(&tmp);
        // a pseudo-random but fixed 900 bp gene
        let mut x: u32 = 12345;
        let gene: Vec<u8> = (0..900)
            .map(|_| {
                x = x.wrapping_mul(1_103_515_245).wrapping_add(12345);
                b"ACGT"[(x >> 16) as usize % 4]
            })
            .collect();
        let mut other = gene.clone();
        other.reverse();
        let filler: Vec<u8> = gene
            .iter()
            .map(|b| match b {
                b'A' => b'C',
                b'C' => b'G',
                b'G' => b'T',
                _ => b'A',
            })
            .collect();
        let mut plasmid = filler.clone();
        plasmid.extend_from_slice(&gene);
        plasmid.extend_from_slice(&filler);
        let write = |name: &str, recs: &[(&str, &[u8])]| {
            let p = tmp.join(name);
            let mut t = Vec::new();
            for (id, s) in recs {
                t.extend_from_slice(format!(">{id}\n").as_bytes());
                t.extend_from_slice(s);
                t.push(b'\n');
            }
            std::fs::write(&p, t).unwrap();
            p
        };
        let mk = |fa: PathBuf, out: &str| {
            let ok = Command::new(&makeblastdb)
                .arg("-in")
                .arg(&fa)
                .args(["-dbtype", "nucl", "-parse_seqids", "-out"])
                .arg(tmp.join("blast").join(out))
                .output()
                .unwrap()
                .status
                .success();
            assert!(ok);
        };
        mk(write("g.fna", &[("NZ_CP2.1", &plasmid)]), "genomes");
        mk(
            write("n.fna", &[("g11", &gene), ("g12", &other)]),
            "groups_nt",
        );

        let hits = lib.search(&blastn, &gene, &tmp.join("work")).unwrap();
        assert_eq!(hits.len(), 1);
        let h = &hits[0];
        assert_eq!(
            (h.accession.as_str(), h.kind.as_str(), h.name.as_str()),
            ("NZ_CP2.1", "plasmid", "pLmN1546")
        );
        assert_eq!((h.start, h.end, h.represents), (901, 1800, 4));
        assert!(h.identity > 99.9 && h.coverage > 99.9, "{h:?}");

        let groups = lib
            .groups_of(&blastn, &[("cadA".into(), gene.clone())], &tmp.join("work"))
            .unwrap();
        assert_eq!(groups.get("cadA"), Some(&vec![11]));
    }

    pub fn tempdir(tag: &str) -> PathBuf {
        let d = std::env::temp_dir().join(format!(
            "sc-library-{tag}-{}-{:?}",
            std::process::id(),
            std::thread::current().id()
        ));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        d
    }
}
