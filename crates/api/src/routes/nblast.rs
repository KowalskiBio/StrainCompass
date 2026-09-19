//! Naming the genes the reference has never seen, through NCBI's public
//! BLAST URL API: the query's novel ORFs are submitted as blastx
//! searches against nr, and the best hit's title becomes the gene's
//! name. This is the only naming resource for truly gained genes - the
//! reference proteome cannot know a gene it does not carry.
//!
//! NCBI asks interactive users of the URL API to poll each search no
//! more often than every few seconds and to keep volumes modest, so the
//! pass is on demand (one region at a time, novel genes only) and its
//! answers are persisted per query: a named gene is never searched
//! twice, and the names surface in the table and exports like the
//! reference-given ones do.

use crate::error::{ApiError, ApiResult};
use crate::state::SharedState;
use axum::extract::{Path, Query, State};
use axum::Json;
use straincompass_types::{GainedRow, IdentifiedOrf, OrfMatch};

/// The BLAST URL API endpoint; overridable so tests can point the whole
/// pass at a stand-in server, the same trade the tool stubs make.
fn blast_url() -> String {
    std::env::var("STRAINCOMPASS_BLAST_URL")
        .unwrap_or_else(|_| "https://blast.ncbi.nlm.nih.gov/Blast.cgi".into())
}

/// A region with more novel genes than this is refused rather than
/// flooding the public service; regions this gene-rich are rare.
const MAX_ORFS: usize = 20;
fn client() -> reqwest::Client {
    reqwest::Client::builder()
        .user_agent(concat!(
            "StrainCompass/",
            env!("CARGO_PKG_VERSION"),
            " (bacterial genome comparison workbench)"
        ))
        .build()
        .expect("reqwest client")
}

/// Submit one blastx search, returning its request id.
async fn submit(c: &reqwest::Client, url: &str, seq: &str) -> ApiResult<String> {
    let resp = c
        .post(url)
        .form(&[
            ("CMD", "Put"),
            ("PROGRAM", "blastx"),
            ("DATABASE", "nr"),
            ("HITLIST_SIZE", "5"),
            ("QUERY", seq),
        ])
        .timeout(std::time::Duration::from_secs(60))
        .send()
        .await
        .map_err(|e| ApiError::Internal(format!("NCBI BLAST could not be reached. ({e})")))?;
    if !resp.status().is_success() {
        return Err(ApiError::Internal(
            "NCBI BLAST refused the search request.".into(),
        ));
    }
    let text = resp
        .text()
        .await
        .map_err(|e| ApiError::Internal(format!("NCBI BLAST answered oddly. ({e})")))?;
    for line in text.lines() {
        if let Some(rid) = line.trim().strip_prefix("RID = ") {
            return Ok(rid.trim().to_string());
        }
    }
    Err(ApiError::Internal(
        "NCBI BLAST did not accept the search.".into(),
    ))
}

/// One poll of one search: `None` while it is still queued, its result
/// text once it is ready.
async fn poll(c: &reqwest::Client, url: &str, rid: &str) -> ApiResult<Option<String>> {
    let resp = c
        .get(url)
        .query(&[
            ("CMD", "Get"),
            ("FORMAT_TYPE", "Text"),
            ("RID", rid),
            ("DESCRIPTIONS", "5"),
            ("ALIGNMENTS", "5"),
        ])
        .timeout(std::time::Duration::from_secs(60))
        .send()
        .await
        .map_err(|e| ApiError::Internal(format!("NCBI BLAST could not be reached. ({e})")))?;
    let text = resp
        .text()
        .await
        .map_err(|e| ApiError::Internal(format!("NCBI BLAST answered oddly. ({e})")))?;
    if text.contains("Status=WAITING")
        || text.contains("Status=QUEUED")
        || text.contains("Status=UNKNOWN")
    {
        if text.contains("Status=UNKNOWN") {
            return Err(ApiError::Internal(
                "NCBI BLAST lost the search before it finished.".into(),
            ));
        }
        return Ok(None);
    }
    Ok(Some(text))
}

/// The best hit of a finished search's text output: title, accession,
/// percent identity, E-value. `None` when nothing significant was
/// found - which is an answer, not an error.
fn parse_best_hit(text: &str) -> Option<(String, String, f64, f64)> {
    // The alignment blocks follow the "Alignments:" header; the first
    // defline is the best hit, and the statistics of its block follow
    // it ("Score = ... Expect = ..." share a line, "Identities = ..."
    // has its own).
    let after = text.split_once("Alignments:")?.1;
    let mut in_block = false;
    let mut title = String::new();
    let mut accession = String::new();
    let mut expect: Option<f64> = None;
    let mut identity: Option<f64> = None;
    for line in after.lines() {
        let line = line.trim();
        if let Some(defline) = line.strip_prefix('>') {
            if in_block {
                // The second defline: the best hit's block is complete.
                break;
            }
            in_block = true;
            let mut parts = defline.splitn(2, char::is_whitespace);
            let raw_id = parts.next().unwrap_or_default();
            // NCBI ids arrive as "db|accession|" or bare accessions.
            let raw_id = raw_id.trim_end_matches('|');
            accession = raw_id.rsplit('|').next().unwrap_or(raw_id).to_string();
            title = parts.next().unwrap_or_default().trim().to_string();
            continue;
        }
        if !in_block {
            continue;
        }
        if expect.is_none() {
            if let Some(pos) = line.find("Expect = ") {
                expect = line[pos + 9..]
                    .split([',', ' '])
                    .next()
                    .and_then(|x| x.trim().parse().ok());
            }
        }
        if identity.is_none() {
            if let Some(v) = line.strip_prefix("Identities = ") {
                identity = v
                    .split('(')
                    .nth(1)
                    .and_then(|p| p.split(')').next())
                    .map(|p| p.trim_end_matches('%'))
                    .and_then(|p| p.parse().ok());
            }
        }
    }
    if title.is_empty() {
        return None;
    }
    Some((
        title,
        accession,
        identity.unwrap_or(0.0),
        expect.unwrap_or(f64::INFINITY),
    ))
}

/// The per-query store of NCBI-given names, one file in the query's
/// directory: region key -> one entry per annotated ORF (index within
/// the row's `orfs`), `null` for "searched, nothing found". Named
/// genes are never searched twice, and the table shows the names
/// without anyone clicking anything again.
fn sidecar_path(qdir: &std::path::Path) -> std::path::PathBuf {
    qdir.join("gained_ncbi_names.json")
}

fn load_sidecar(qdir: &std::path::Path) -> serde_json::Map<String, serde_json::Value> {
    std::fs::read(sidecar_path(qdir))
        .ok()
        .and_then(|b| serde_json::from_slice(&b).ok())
        .unwrap_or_default()
}

fn save_sidecar(qdir: &std::path::Path, v: &serde_json::Map<String, serde_json::Value>) {
    if let Ok(b) = serde_json::to_vec_pretty(v) {
        let _ = std::fs::write(sidecar_path(qdir), b);
    }
}

/// Fill each row's ORF `ncbi` field from the sidecar, so the table and
/// the exports carry the NCBI-given names the same way they carry the
/// reference-given ones. Called wherever gained rows are read.
pub fn merge_ncbi_names(qdir: &std::path::Path, rows: &mut [GainedRow]) {
    if rows.is_empty() {
        return;
    }
    let sidecar = load_sidecar(qdir);
    for row in rows.iter_mut() {
        let Some(entries) = sidecar.get(&format!("{}:{}-{}", row.qry_seqid, row.start, row.end))
        else {
            continue;
        };
        let Some(list) = entries.as_array() else {
            continue;
        };
        for e in list {
            let idx = e.get("orf").and_then(|v| v.as_u64()).unwrap_or(u64::MAX) as usize;
            if let Some(orf) = row.orfs.get_mut(idx) {
                orf.ncbi = e
                    .get("match")
                    .and_then(|m| serde_json::from_value::<OrfMatch>(m.clone()).ok());
            }
        }
    }
}

/// The query-agnostic shape of one sidecar entry. `source` tells the
/// on-demand nr pass which answers are its own: a null from SwissProt
/// still leaves the gene worth an nr search, a null from nr is final.
fn entry(orf: usize, m: &Option<OrfMatch>, source: &str) -> serde_json::Value {
    serde_json::json!({
        "orf": orf,
        "match": m,
        "source": source,
    })
}

// ---------------------------------------------------------------------
// The automatic naming pass
//
// A finished comparison names its novel genes right away, against a
// local copy of NCBI's curated SwissProt database: searching the public
// BLAST service for a whole run's worth of genes (thousands) would
// abuse it and take days, while the local search takes minutes and
// names the classes that matter (mobilization, phage, resistance).
// Genes SwissProt misses keep the on-demand nr search from the region
// card.
// ---------------------------------------------------------------------

/// Where the curated database lives: `STRAINCOMPASS_SWISSPROT_DB`, or
/// the conventional `blastdb/swissprot` beside the data directory. A
/// missing database disables the pass quietly - deployments made
/// before it existed lose nothing.
fn swissprot_db(state: &SharedState) -> Option<std::path::PathBuf> {
    let p = std::env::var("STRAINCOMPASS_SWISSPROT_DB")
        .map(std::path::PathBuf::from)
        .ok()
        .or_else(|| {
            state
                .data_dir
                .parent()
                .map(|p| p.join("blastdb").join("swissprot"))
        })?;
    if p.with_extension("pin").is_file() {
        Some(p)
    } else {
        None
    }
}

/// The pass's progress, persisted in the run directory so the UI can
/// pick it up across reloads and restarts.
#[derive(serde::Serialize, serde::Deserialize, Clone)]
pub struct NcbiStatus {
    pub state: String,
    pub named: usize,
    pub total: usize,
    /// Unix seconds of the last update; a "running" status older than
    /// a quarter hour means the pass died (server restart, say) and
    /// must not spin the badge forever.
    pub updated: i64,
}

fn write_status(state: &SharedState, project_id: i64, run_id: i64, s: &NcbiStatus) {
    if let Ok(b) = serde_json::to_vec(s) {
        let _ = std::fs::write(
            state.run_dir(project_id, run_id).join("ncbi_status.json"),
            b,
        );
    }
}

/// GET /runs/{id}/ncbi_status : how the run's automatic naming is
/// getting on.
pub async fn gained_ncbi_status(
    State(state): State<SharedState>,
    Path(run_id): Path<i64>,
) -> ApiResult<Json<NcbiStatus>> {
    let (project_id, _, _) = crate::jobs::run_meta(&state, run_id)?;
    let file = state.run_dir(project_id, run_id).join("ncbi_status.json");
    let mut s: NcbiStatus = std::fs::read(&file)
        .ok()
        .and_then(|b| serde_json::from_slice(&b).ok())
        .unwrap_or(NcbiStatus {
            state: "idle".into(),
            named: 0,
            total: 0,
            updated: 0,
        });
    if s.state == "running" && s.updated > 0 {
        let now = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_secs() as i64)
            .unwrap_or(0);
        if now - s.updated > 900 {
            s.state = "interrupted".into();
        }
    }
    Ok(Json(s))
}

/// One sidecar entry, parsed: the ORF index it belongs to, the name it
/// carries (if any), and which pass gave it.
#[derive(Clone)]
struct SidecarEntry {
    orf: usize,
    m: Option<OrfMatch>,
    source: String,
}

fn parse_entries(v: Option<&serde_json::Value>) -> Vec<SidecarEntry> {
    v.and_then(|v| v.as_array())
        .map(|list| {
            list.iter()
                .filter_map(|e| {
                    Some(SidecarEntry {
                        orf: e.get("orf")?.as_u64()? as usize,
                        m: e.get("match")
                            .and_then(|m| serde_json::from_value::<OrfMatch>(m.clone()).ok()),
                        source: e
                            .get("source")
                            .and_then(|s| s.as_str())
                            .unwrap_or("nr")
                            .to_string(),
                    })
                })
                .collect()
        })
        .unwrap_or_default()
}

/// The automatic pass itself, spawned when a comparison succeeds.
/// One query at a time: load its gained rows, find the genes with no
/// name of any kind, answer from the project-wide cache where
/// possible, search the rest against SwissProt, persist everything,
/// and heartbeat the status for the badge.
async fn auto_pass(
    state: &SharedState,
    project_id: i64,
    run_id: i64,
    query_ids: &[i64],
) -> ApiResult<()> {
    let Some(db) = swissprot_db(state) else {
        return Ok(());
    };
    let run_dir = state.run_dir(project_id, run_id);

    // Count the work first, so the badge can say what it is waiting
    // for.
    let mut total = 0usize;
    let mut per_query: Vec<Vec<(usize, usize)>> = Vec::new(); // (row, orf) indices
    let mut rows_all: Vec<Vec<straincompass_types::GainedRow>> = Vec::new();
    for qid in query_ids {
        let res = crate::jobs::load_query_result(state, project_id, run_id, *qid)?;
        let rows = res.gained.unwrap_or_default();
        // A gene with any stored answer - a name from any source, an
        // nr "nothing found", or a search still queued - is not this
        // pass's to search.
        let sidecar = load_sidecar(&run_dir.join("queries").join(qid.to_string()));
        let mut novel = Vec::new();
        for (ri, r) in rows.iter().enumerate() {
            let entries =
                parse_entries(sidecar.get(&format!("{}:{}-{}", r.qry_seqid, r.start, r.end)));
            for (oi, o) in r.orfs.iter().enumerate() {
                if o.best.is_none() && !entries.iter().any(|e| e.orf == oi) {
                    novel.push((ri, oi));
                }
            }
        }
        total += novel.len();
        per_query.push(novel);
        rows_all.push(rows);
    }
    write_status(
        state,
        project_id,
        run_id,
        &NcbiStatus {
            state: "running".into(),
            named: 0,
            total,
            updated: now_secs(),
        },
    );
    if total == 0 {
        write_status(
            state,
            project_id,
            run_id,
            &NcbiStatus {
                state: "done".into(),
                named: 0,
                total: 0,
                updated: now_secs(),
            },
        );
        return Ok(());
    }

    // The project-wide cache: a gene's sequence answers the same
    // wherever it appears, so re-runs and neighbouring strains cost
    // nothing. Keyed by the DNA's digest.
    let cache_path = state.project_dir(project_id).join("ncbi_cache.json");
    let mut cache: std::collections::HashMap<String, Option<OrfMatch>> = std::fs::read(&cache_path)
        .ok()
        .and_then(|b| serde_json::from_slice(&b).ok())
        .unwrap_or_default();

    let mut named_total = 0usize;
    for (qi, qid) in query_ids.iter().enumerate() {
        let qdir = run_dir.join("queries").join(qid.to_string());
        let novel = &per_query[qi];
        if novel.is_empty() {
            continue;
        }

        // The genes' sequences, from the query fasta.
        let jobs: Vec<((usize, usize), String)> = {
            let qry_fa = qdir.join("query.fa");
            let ids: Vec<(usize, usize)> = novel.clone();
            let rows_c = rows_all[qi].clone();
            let seqs = tokio::task::spawn_blocking(move || {
                let records = straincompass_engine::fasta::parse_fasta(&qry_fa)?;
                let mut out = Vec::new();
                for (ri, oi) in ids {
                    out.push(rows_all_snapshot(&rows_c, ri, oi, &records)?);
                }
                Ok::<_, ApiError>(out)
            })
            .await
            .map_err(|e| ApiError::Internal(format!("Reading the genes crashed. ({e})")))??;
            novel.iter().cloned().zip(seqs).collect()
        };

        // Cached answers skip the search entirely.
        let mut to_search: Vec<((usize, usize), String)> = Vec::new();
        let mut results: Vec<((usize, usize), Option<OrfMatch>)> = Vec::new();
        for (p, seq) in jobs {
            let key = digest(&seq);
            if let Some(m) = cache.get(&key) {
                results.push((p, m.clone()));
            } else {
                to_search.push((p, seq));
            }
        }

        if !to_search.is_empty() {
            let db = db.clone();
            let searches = to_search.clone();
            let found = tokio::task::spawn_blocking(move || search_swissprot(&db, &searches))
                .await
                .map_err(|e| ApiError::Internal(format!("The naming search crashed. ({e})")))??;
            for ((p, seq), m) in to_search.into_iter().zip(found.iter()) {
                cache.insert(digest(&seq), m.clone());
                results.push((p, m.clone()));
            }
            let _ = std::fs::write(&cache_path, serde_json::to_vec(&cache).unwrap_or_default());
        }

        // Persist into the same sidecar the on-demand pass uses.
        let mut sidecar = load_sidecar(&qdir);
        for ((ri, oi), m) in results {
            if m.is_some() {
                named_total += 1;
            }
            let row = &rows_all[qi][ri];
            let key = format!("{}:{}-{}", row.qry_seqid, row.start, row.end);
            let mut entries = parse_entries(sidecar.get(&key));
            entries.retain(|e| e.orf != oi);
            entries.push(SidecarEntry {
                orf: oi,
                m,
                source: "swissprot".into(),
            });
            let list: Vec<serde_json::Value> = entries
                .iter()
                .map(|e| entry(e.orf, &e.m, &e.source))
                .collect();
            sidecar.insert(key, serde_json::Value::Array(list));
        }
        save_sidecar(&qdir, &sidecar);

        write_status(
            state,
            project_id,
            run_id,
            &NcbiStatus {
                state: "running".into(),
                named: named_total,
                total,
                updated: now_secs(),
            },
        );
    }

    write_status(
        state,
        project_id,
        run_id,
        &NcbiStatus {
            state: "done".into(),
            named: named_total,
            total,
            updated: now_secs(),
        },
    );
    Ok(())
}

fn now_secs() -> i64 {
    std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_secs() as i64)
        .unwrap_or(0)
}

fn rows_all_snapshot(
    rows: &[straincompass_types::GainedRow],
    ri: usize,
    oi: usize,
    records: &[straincompass_engine::fasta::FastaRecord],
) -> ApiResult<String> {
    let row = &rows[ri];
    let o = &row.orfs[oi];
    let rec = records
        .iter()
        .find(|r| r.id == row.qry_seqid)
        .ok_or_else(|| {
            ApiError::BadRequest("The query contig is not in this query's fasta file.".into())
        })?;
    Ok(
        String::from_utf8_lossy(&straincompass_engine::fasta::subseq(
            rec,
            o.start,
            o.end,
            o.strand < 0,
        ))
        .to_ascii_uppercase(),
    )
}

fn digest(s: &str) -> String {
    use sha2::Digest as _;
    let d = sha2::Sha256::digest(s.as_bytes());
    d.iter().map(|b| format!("{b:02x}")).collect()
}

/// One blastx of every uncached gene against the curated database,
/// tabular, best hit per gene by bitscore.
fn search_swissprot(
    db: &std::path::Path,
    searches: &[((usize, usize), String)],
) -> ApiResult<Vec<Option<OrfMatch>>> {
    use std::io::Write as _;
    let tools = straincompass_engine::tools::ToolPaths::discover()
        .map_err(|e| ApiError::Internal(e.to_string()))?;
    let tmp = std::env::temp_dir().join(format!(
        "straincompass-ncbi-{}-{:x}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_nanos())
            .unwrap_or(0)
    ));
    std::fs::create_dir_all(&tmp)?;
    let fa = tmp.join("orfs.fa");
    {
        let mut f = std::fs::File::create(&fa)?;
        for (i, ((_, _), seq)) in searches.iter().enumerate() {
            writeln!(f, ">g{i}\n{seq}")?;
        }
    }
    let out = std::process::Command::new(&tools.blastx)
        .arg("-query")
        .arg(&fa)
        .arg("-db")
        .arg(db)
        .arg("-evalue")
        .arg("1e-5")
        .arg("-max_hsps")
        .arg("1")
        .arg("-num_threads")
        .arg("4")
        .arg("-outfmt")
        .arg("6 qseqid pident qcovhsp evalue bitscore stitle")
        .output()
        .map_err(|e| ApiError::Internal(format!("blastx could not be run. ({e})")))?;
    let _ = std::fs::remove_dir_all(&tmp);
    if !out.status.success() {
        return Err(ApiError::Internal(format!(
            "blastx failed: {}",
            String::from_utf8_lossy(&out.stderr)
        )));
    }

    // Best hit per gene, by bitscore.
    let mut best: std::collections::HashMap<String, OrfMatch> = std::collections::HashMap::new();
    for line in String::from_utf8_lossy(&out.stdout).lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 6 {
            continue;
        }
        let bitscore: f64 = match f[4].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let better = best
            .get(f[0])
            .map(|m: &OrfMatch| bitscore > m.bitscore)
            .unwrap_or(true);
        if !better {
            continue;
        }
        // "sp|P0A3H9|RELX_LISMO Relaxase/mobilization protein n..." -
        // the accession sits in the id, the name after it.
        let mut title_parts = f[5].splitn(2, char::is_whitespace);
        let raw_id = title_parts.next().unwrap_or_default();
        let accession = raw_id
            .trim_end_matches('|')
            .split('|')
            .nth(1)
            .unwrap_or(raw_id)
            .to_string();
        let label = title_parts.next().unwrap_or_default().trim().to_string();
        best.insert(
            f[0].to_string(),
            OrfMatch {
                locus_tag: accession.clone(),
                protein_id: accession,
                label,
                identity: f[1].parse().unwrap_or(0.0),
                coverage: f[2].parse().unwrap_or(0.0),
                evalue: f[3].parse().unwrap_or(f64::INFINITY),
                bitscore,
            },
        );
    }
    Ok((0..searches.len())
        .map(|i| best.get(&format!("g{i}")).cloned())
        .collect())
}

/// Called when a comparison succeeds: name its novel genes in the
/// background. The comparison's own speed is untouched - the pass
/// works from the finished result files.
pub fn spawn_auto_pass(state: SharedState, project_id: i64, run_id: i64, query_ids: Vec<i64>) {
    tokio::spawn(async move {
        if let Err(e) = auto_pass(&state, project_id, run_id, &query_ids).await {
            tracing::error!("run {run_id}: the automatic naming failed: {e}");
            write_status(
                &state,
                project_id,
                run_id,
                &NcbiStatus {
                    state: "failed".into(),
                    named: 0,
                    total: 0,
                    updated: now_secs(),
                },
            );
        }
    });
}
/// The state of one region's NCBI naming, returned to the client: the
/// ORFs with whatever names are known so far, and how many searches are
/// still queued at NCBI (the client polls again until it reaches
/// zero).
#[derive(serde::Serialize)]
pub struct GainedNcbiStatus {
    pub orfs: Vec<IdentifiedOrf>,
    pub pending: u32,
}

/// GET /runs/{id}/gained/annotate_ncbi?query_id&seqid&start&end : name
/// one region's novel genes through NCBI BLAST.
///
/// The service's queue is minutes, not seconds, so this is a
/// submit-and-poll dance the client drives: the first call submits the
/// region's unnamed genes and returns immediately; each later call
/// polls the outstanding searches once and returns what has finished.
/// The answers are persisted as they arrive, so closing the card
/// mid-search loses nothing - the next call picks the searches up, and
/// the table shows the names from then on.
pub async fn gained_annotate_ncbi(
    State(state): State<SharedState>,
    Path(run_id): Path<i64>,
    Query(q): Query<super::results::GainedVerifyQuery>,
) -> ApiResult<Json<GainedNcbiStatus>> {
    let (project_id, run_id, qid, _, _, _, row) =
        super::results::gained_region_ctx(&state, run_id, &q)?;
    let qdir = state
        .run_dir(project_id, run_id)
        .join("queries")
        .join(qid.to_string());

    let mut sidecar = load_sidecar(&qdir);
    let region_key = format!("{}:{}-{}", row.qry_seqid, row.start, row.end);
    // The stored entries, by what they mean for the nr search: a gene
    // with a name of any source (reference, SwissProt, nr) needs
    // nothing; an nr answer - a name or "nothing found" - is final for
    // this pass; a SwissProt miss still deserves the nr search; a "rid"
    // is an nr search still queued at NCBI.
    let mut entries = parse_entries(sidecar.get(&region_key));
    entries.retain(|e| e.orf < row.orfs.len());
    let mut named: std::collections::HashMap<usize, OrfMatch> = entries
        .iter()
        .filter_map(|e| e.m.clone().map(|m| (e.orf, m)))
        .collect();
    let mut nr_done: std::collections::HashSet<usize> = entries
        .iter()
        .filter(|e| e.source == "nr" && e.m.is_none())
        .map(|e| e.orf)
        .collect();
    let mut pending: Vec<(usize, String)> = Vec::new();
    for e in sidecar
        .get(&region_key)
        .and_then(|v| v.as_array())
        .into_iter()
        .flatten()
    {
        if let Some(rid) = e.get("rid").and_then(|v| v.as_str()) {
            if let Some(idx) = e.get("orf").and_then(|v| v.as_u64()) {
                pending.push((idx as usize, rid.to_string()));
            }
        }
    }

    // The genes worth an nr submission: novel, unnamed by any source,
    // never searched by nr, not queued.
    let fresh: Vec<usize> = (0..row.orfs.len())
        .filter(|i| {
            row.orfs[*i].best.is_none()
                && !named.contains_key(i)
                && !nr_done.contains(i)
                && !pending.iter().any(|(j, _)| j == i)
        })
        .collect();
    if entries.len() + fresh.len() > MAX_ORFS {
        return Err(ApiError::BadRequest(format!(
            "This region has {} genes to name, more than the {} the NCBI service should be asked about at once. Search the region at NCBI BLAST from the region card instead.",
            entries.len() + fresh.len(),
            MAX_ORFS
        )));
    }

    if !fresh.is_empty() {
        let qry_fa = qdir.join("query.fa");
        let novel_idx = fresh.clone();
        let row_c = row.clone();
        let seqs = tokio::task::spawn_blocking(move || {
            let records = straincompass_engine::fasta::parse_fasta(&qry_fa)?;
            let rec = records
                .iter()
                .find(|r| r.id == row_c.qry_seqid)
                .ok_or_else(|| {
                    ApiError::BadRequest(
                        "The query contig is not in this query's fasta file.".into(),
                    )
                })?;
            let out: Vec<String> = novel_idx
                .iter()
                .map(|i| {
                    let o = &row_c.orfs[*i];
                    String::from_utf8_lossy(&straincompass_engine::fasta::subseq(
                        rec,
                        o.start,
                        o.end,
                        o.strand < 0,
                    ))
                    .to_ascii_uppercase()
                })
                .collect();
            Ok::<_, ApiError>(out)
        })
        .await
        .map_err(|e| ApiError::Internal(format!("Reading the genes crashed. ({e})")))??;

        let url = blast_url();
        let c = client();
        for (i, seq) in fresh.iter().zip(seqs.iter()) {
            let rid = submit(&c, &url, seq).await?;
            pending.push((*i, rid));
        }
    }

    // One poll round of every outstanding search; finished ones become
    // answers. A search NCBI has lost is dropped, so a later call can
    // submit it again rather than waiting forever.
    if !pending.is_empty() {
        let url = blast_url();
        let c = client();
        let mut still: Vec<(usize, String)> = Vec::new();
        for (i, rid) in pending {
            match poll(&c, &url, &rid).await? {
                Some(text) => {
                    let m = parse_best_hit(&text).map(|(label, protein_id, identity, evalue)| {
                        OrfMatch {
                            locus_tag: protein_id.clone(),
                            protein_id,
                            label,
                            identity,
                            coverage: 0.0,
                            evalue,
                            bitscore: 0.0,
                        }
                    });
                    if let Some(m) = m {
                        named.insert(i, m);
                    } else {
                        nr_done.insert(i);
                    }
                }
                None => still.push((i, rid)),
            }
        }
        pending = still;
    }

    // Persist before answering: submissions must not be lost, and the
    // names must survive the session. Prior entries of other sources
    // stay - the automatic pass's names are as good as these.
    let list: Vec<serde_json::Value> = entries
        .iter()
        .filter(|e| !named.contains_key(&e.orf) && !nr_done.contains(&e.orf))
        .map(|e| entry(e.orf, &e.m, &e.source))
        .chain(named.iter().map(|(i, m)| entry(*i, &Some(m.clone()), "nr")))
        .chain(nr_done.iter().map(|i| entry(*i, &None, "nr")))
        .chain(
            pending
                .iter()
                .map(|(i, rid)| serde_json::json!({ "orf": i, "rid": rid })),
        )
        .collect();
    sidecar.insert(region_key, serde_json::Value::Array(list));
    save_sidecar(&qdir, &sidecar);

    let orfs: Vec<IdentifiedOrf> = row
        .orfs
        .iter()
        .enumerate()
        .map(|(i, o)| IdentifiedOrf {
            start: o.start,
            end: o.end,
            strand: o.strand,
            best: named.get(&i).cloned(),
            seq: String::new(),
        })
        .collect();
    Ok(Json(GainedNcbiStatus {
        orfs,
        pending: pending.len() as u32,
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    const TEXT: &str = "\
BLASTX 2.16.0+\n\
\n\
Query= o3\n\
Length=729\n\
\n\
Sequences producing significant alignments:\n\
Score E Sequences producing significant alignments\n\
\n\
Alignments:\n\
>gb|WP_097529575.1| MobV family relaxase, partial [Listeria monocytogenes]\n\
Length=285\n\n\
 Score = 521 bits (1342),  Expect = 1e-170, Method: Compositional matrix adjust.\n\
 Identities = 253/253 (100%), Gaps = 0/253 (0%)\n\
 Frame = +1\n\
\n\
>gb|WP_012680999.1| chromosomal replication initiator protein DnaA\n\
Length=452\n\n\
 Score = 32 bits (70),  Expect = 1e-5, Method: Compositional matrix adjust.\n\
 Identities = 15/52 (28%), Gaps = 0/52 (0%)\n\
 Frame = +3\n\
";

    #[test]
    fn the_best_hit_is_parsed_with_its_statistics() {
        let (label, accession, identity, evalue) =
            parse_best_hit(TEXT).expect("the first hit must be read");
        assert_eq!(
            label,
            "MobV family relaxase, partial [Listeria monocytogenes]"
        );
        assert_eq!(accession, "WP_097529575.1");
        assert_eq!(identity, 100.0);
        assert_eq!(evalue, 1e-170);
    }

    #[test]
    fn no_hits_is_not_an_error() {
        let text = "BLASTX 2.16.0+\n\n\nNo significant similarity found.\n";
        assert!(parse_best_hit(text).is_none());
    }
}

#[cfg(test)]
mod endpoint_tests {
    use super::*;
    use axum::extract::{Path, Query, State};

    /// A stand-in NCBI BLAST service: every POST is a submission and
    /// gets a fixed RID, every GET is a poll and gets a finished search
    /// with one MobV-shaped hit. Speaking just enough HTTP for reqwest.
    async fn stub_ncbi(listener: tokio::net::TcpListener) {
        use tokio::io::{AsyncReadExt, AsyncWriteExt};
        loop {
            let (mut sock, _) = match listener.accept().await {
                Ok(s) => s,
                Err(_) => return,
            };
            tokio::spawn(async move {
                let mut buf = Vec::new();
                let mut tmp = [0u8; 4096];
                // Read headers (and whatever body arrives with them);
                // good enough for the small requests the pass makes.
                loop {
                    match sock.read(&mut tmp).await {
                        Ok(0) | Err(_) => return,
                        Ok(n) => {
                            buf.extend_from_slice(&tmp[..n]);
                            if buf.windows(4).any(|w| w == b"\r\n\r\n") {
                                break;
                            }
                        }
                    }
                }
                let head = String::from_utf8_lossy(&buf).to_string();
                let is_get = head.starts_with("GET");
                let body = if is_get {
                    "HTTP/1.1 200 OK\r\nContent-Type: text/plain\r\n\r\n\
                     BLASTX 2.16.0+\n\n\nAlignments:\n\
                     >gb|WP_097529575.1| MobV family relaxase, partial [Listeria monocytogenes]\n\
                     Length=285\n\n\
                     Score = 521 bits (1342),  Expect = 1e-170\n\
                     Identities = 253/253 (100%), Gaps = 0/253 (0%)\n\n"
                        .to_string()
                } else {
                    "HTTP/1.1 200 OK\r\nContent-Type: text/plain\r\n\r\nRID = STUBRID01\nRTOE = 1\n"
                        .to_string()
                };
                let _ = sock.write_all(body.as_bytes()).await;
                let _ = sock.shutdown().await;
            });
        }
    }

    #[tokio::test]
    async fn novel_genes_are_named_and_persisted() {
        let listener = tokio::net::TcpListener::bind("127.0.0.1:0").await.unwrap();
        let addr = listener.local_addr().unwrap();
        tokio::spawn(stub_ncbi(listener));
        std::env::set_var(
            "STRAINCOMPASS_BLAST_URL",
            format!("http://{addr}/Blast.cgi"),
        );

        let (state, dir) = crate::routes::results::tests::seeded_state();
        crate::routes::results::tests::seed_gained_on_real_contig(&dir);

        let orfs = gained_annotate_ncbi(
            State(state.clone()),
            Path(1),
            Query(crate::routes::results::GainedVerifyQuery {
                query_id: None,
                seqid: "q1".into(),
                start: 1,
                end: 6,
            }),
        )
        .await
        .unwrap()
        .0;
        // The seeded ORF has no reference name, so the stub service gets
        // to name it. The stub answers polls immediately, so one call
        // submits and collects.
        assert_eq!(orfs.pending, 0, "nothing left queued at the stub");
        let orfs = &orfs.orfs;
        assert_eq!(orfs.len(), 1);
        let m = orfs[0].best.as_ref().expect("named by the stub service");
        assert_eq!(
            m.label,
            "MobV family relaxase, partial [Listeria monocytogenes]"
        );
        assert_eq!(m.protein_id, "WP_097529575.1");
        assert_eq!(m.identity, 100.0);

        // The answer is persisted, so the table read merges it and the
        // second call costs no submission at all.
        let sidecar = std::fs::read_to_string(
            dir.join("projects/1/runs/1/queries/10/gained_ncbi_names.json"),
        )
        .unwrap();
        assert!(sidecar.contains("MobV family relaxase"));

        let orfs2 = gained_annotate_ncbi(
            State(state),
            Path(1),
            Query(crate::routes::results::GainedVerifyQuery {
                query_id: None,
                seqid: "q1".into(),
                start: 1,
                end: 6,
            }),
        )
        .await
        .unwrap()
        .0;
        // No new submissions were made (the stub RID appears once), the
        // stored answer is returned as-is.
        assert_eq!(orfs2.pending, 0);
        assert_eq!(orfs2.orfs[0].best.as_ref().unwrap().label, m.label);

        std::env::remove_var("STRAINCOMPASS_BLAST_URL");
        let _ = std::fs::remove_dir_all(&dir);
    }
}

#[cfg(test)]
mod auto_pass_tests {
    use super::*;

    /// The automatic pass, end to end against a stand-in curated
    /// database: the novel gene is searched once, named, persisted in
    /// the sidecar, and the status the badge polls says so.
    #[tokio::test]
    async fn the_auto_pass_names_novel_genes_and_persists() {
        let _env = crate::routes::results::tests::TOOLS_ENV.lock().await;
        let (state, dir) = crate::routes::results::tests::seeded_state();
        crate::routes::results::tests::seed_gained_on_real_contig(&dir);
        crate::routes::results::tests::stub_blast_tools(&dir);
        // Its blastx writes to -out; this pass reads stdout, so give
        // the stub a stdout-speaking replacement.
        {
            use std::os::unix::fs::PermissionsExt;
            let blastx = dir.join("stubbin/blastx");
            std::fs::write(
                &blastx,
                "#!/bin/sh\nprintf 'g0\\t99.0\\t100.0\\t1e-180\\t520\\tsp|P0A3H9|RELX_LISMO Relaxase/mobilization protein n\\n'\n",
            )
            .unwrap();
            std::fs::set_permissions(&blastx, std::fs::Permissions::from_mode(0o755)).unwrap();
        }
        // A stand-in database: the pass only checks the index exists.
        let db = dir.join("blastdb/swissprot");
        std::fs::create_dir_all(&db).unwrap();
        std::fs::write(db.with_extension("pin"), b"stub").unwrap();
        std::env::set_var("STRAINCOMPASS_SWISSPROT_DB", &db);

        auto_pass(&state, 1, 1, &[10]).await.unwrap();

        // The status the badge polls: done, the one novel gene named.
        let st = gained_ncbi_status(State(state.clone()), Path(1))
            .await
            .unwrap()
            .0;
        assert_eq!((st.state.as_str(), st.named, st.total), ("done", 1, 1));

        // The answer rides the sidecar, and the table read merges it
        // like any other NCBI-given name.
        let qdir = dir.join("projects/1/runs/1/queries/10");
        let sidecar = std::fs::read_to_string(qdir.join("gained_ncbi_names.json")).unwrap();
        assert!(sidecar.contains("swissprot"));
        assert!(sidecar.contains("Relaxase/mobilization protein n"));
        let mut rows = crate::jobs::load_query_result(&state, 1, 1, 10)
            .unwrap()
            .gained
            .unwrap();
        merge_ncbi_names(&qdir, &mut rows);
        assert_eq!(rows[0].orfs[0].ncbi.as_ref().unwrap().protein_id, "P0A3H9");

        // The project-wide cache holds the answer for every future run.
        assert!(dir.join("projects/1/ncbi_cache.json").is_file());

        std::env::remove_var("STRAINCOMPASS_SWISSPROT_DB");
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// Without the curated database the pass is a quiet no-op: nothing
    /// written, nothing thrown.
    #[tokio::test]
    async fn without_the_database_the_pass_is_a_noop() {
        let _env = crate::routes::results::tests::TOOLS_ENV.lock().await;
        let (state, dir) = crate::routes::results::tests::seeded_state();
        crate::routes::results::tests::seed_gained_on_real_contig(&dir);
        std::env::remove_var("STRAINCOMPASS_SWISSPROT_DB");

        auto_pass(&state, 1, 1, &[10]).await.unwrap();
        assert!(!dir.join("projects/1/runs/1/ncbi_status.json").is_file());

        let _ = std::fs::remove_dir_all(&dir);
    }
}
