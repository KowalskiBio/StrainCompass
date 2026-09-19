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
/// Poll interval and budget: NCBI's queue for a short blastx is tens of
/// seconds, four minutes covers the slow days. The interval bends to
/// the environment so the stub-server test does not wait for it.
const MAX_POLLS: u32 = 20;

fn poll_secs() -> u64 {
    std::env::var("STRAINCOMPASS_TEST_POLL_SECS")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(12)
}

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

/// The query-agnostic shape of one sidecar entry.
fn entry(orf: usize, m: &Option<OrfMatch>) -> serde_json::Value {
    serde_json::json!({
        "orf": orf,
        "match": m,
    })
}

/// GET /runs/{id}/gained/annotate_ncbi?query_id&seqid&start&end : name
/// one region's novel genes through NCBI BLAST, persist the answer, and
/// return the region's ORFs with whatever names are now known.
#[allow(clippy::too_many_arguments)]
pub async fn gained_annotate_ncbi(
    State(state): State<SharedState>,
    Path(run_id): Path<i64>,
    Query(q): Query<super::results::GainedVerifyQuery>,
) -> ApiResult<Json<Vec<IdentifiedOrf>>> {
    let (project_id, run_id, qid, _, _, _, row) =
        super::results::gained_region_ctx(&state, run_id, &q)?;
    let qdir = state
        .run_dir(project_id, run_id)
        .join("queries")
        .join(qid.to_string());

    // The novel ORFs: no reference protein and no stored NCBI answer.
    let mut sidecar = load_sidecar(&qdir);
    let region_key = format!("{}:{}-{}", row.qry_seqid, row.start, row.end);
    let mut stored: std::collections::HashMap<usize, Option<OrfMatch>> = sidecar
        .get(&region_key)
        .and_then(|v| v.as_array())
        .map(|list| {
            list.iter()
                .filter_map(|e| {
                    let idx = e.get("orf")?.as_u64()? as usize;
                    let m = e
                        .get("match")
                        .and_then(|m| serde_json::from_value::<OrfMatch>(m.clone()).ok());
                    Some((idx, m))
                })
                .collect()
        })
        .unwrap_or_default();

    let novel: Vec<usize> = (0..row.orfs.len())
        .filter(|i| row.orfs[*i].best.is_none() && !stored.contains_key(i))
        .collect();
    if novel.len() > MAX_ORFS {
        return Err(ApiError::BadRequest(format!(
            "This region has {} genes without names, more than the {} the NCBI service should be asked about at once. Search the region at NCBI BLAST from the region card instead.",
            novel.len(),
            MAX_ORFS
        )));
    }

    if !novel.is_empty() {
        // The sequences come from the query fasta; the region's row has
        // the coordinates.
        let qry_fa = state
            .run_dir(project_id, run_id)
            .join("queries")
            .join(qid.to_string())
            .join("query.fa");
        let novel_idx = novel.clone();
        let row_c = row.clone();
        let seqs = tokio::task::spawn_blocking(move || {
            let records = straincompass_engine::fasta::parse_fasta(&qry_fa)?;
            let rec = records
                .iter()
                .find(|r| r.id == row_c.qry_seqid)
                .ok_or_else(|| {
                    crate::error::ApiError::BadRequest(
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
            Ok::<_, crate::error::ApiError>(out)
        })
        .await
        .map_err(|e| ApiError::Internal(format!("Reading the genes crashed. ({e})")))??;

        let url = blast_url();
        let c = client();
        let mut rids: Vec<(usize, String)> = Vec::new();
        for (i, seq) in novel.iter().zip(seqs.iter()) {
            rids.push((*i, submit(&c, &url, seq).await?));
        }

        let mut results: Vec<(usize, Option<OrfMatch>)> = Vec::new();
        let mut pending = rids;
        for attempt in 0..MAX_POLLS {
            if pending.is_empty() {
                break;
            }
            tokio::time::sleep(std::time::Duration::from_secs(poll_secs())).await;
            let mut still = Vec::new();
            for (i, rid) in pending {
                match poll(&c, &url, &rid).await? {
                    Some(text) => {
                        let m =
                            parse_best_hit(&text).map(|(label, protein_id, identity, evalue)| {
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
                        results.push((i, m));
                    }
                    None => still.push((i, rid)),
                }
            }
            pending = still;
            if attempt == MAX_POLLS - 1 && !pending.is_empty() {
                return Err(ApiError::Internal(
                    "NCBI BLAST did not finish the search in time. Try again in a moment.".into(),
                ));
            }
        }
        for (i, m) in results {
            stored.insert(i, m);
        }
        // Persist before answering: the names must survive the session.
        let list: Vec<serde_json::Value> = stored.iter().map(|(i, m)| entry(*i, m)).collect();
        sidecar.insert(region_key.clone(), serde_json::Value::Array(list));
        save_sidecar(&qdir, &sidecar);
    }

    let orfs: Vec<IdentifiedOrf> = row
        .orfs
        .iter()
        .enumerate()
        .map(|(i, o)| IdentifiedOrf {
            start: o.start,
            end: o.end,
            strand: o.strand,
            best: stored.get(&i).cloned().flatten(),
            seq: String::new(),
        })
        .collect();
    Ok(Json(orfs))
}

/// Tests for the parser against a trimmed but faithful slice of NCBI's
/// text output.
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
        std::env::set_var("STRAINCOMPASS_TEST_POLL_SECS", "0");

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
        // to name it.
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
        assert_eq!(orfs2[0].best.as_ref().unwrap().label, m.label);

        std::env::remove_var("STRAINCOMPASS_BLAST_URL");
        std::env::remove_var("STRAINCOMPASS_TEST_POLL_SECS");
        let _ = std::fs::remove_dir_all(&dir);
    }
}
