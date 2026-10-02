//! Where a panel gene usually occurs: an NCBI BLAST search of the gene's
//! sequence, its full-length matches sorted into plasmid, chromosome and
//! unplaced-contig records. It answers what the run cannot when no strain
//! carries the gene (bcrA on pLM80-type plasmids), and hands the user a
//! complete plasmid to run the whole-element comparison with.
//!
//! Only the panel gene's own sequence is sent to NCBI, never a strain's
//! genome. The search is scoped to the project's genus (from the
//! project's organism, else the reference FASTA header) and widened to
//! all bacteria when the genus has too few matches; the answer always
//! states its scope.
//!
//! Like the novel-gene naming, nothing runs in the background: the first
//! request submits the search and stores its request id, and each later
//! request polls it once - at most once a minute, as NCBI asks - until
//! the answer is stored for good.

use crate::error::{ApiError, ApiResult};
use crate::jobs;
use crate::routes::nblast::{blast_url, client};
use crate::state::SharedState;
use axum::extract::{Path, Query, State};
use axum::Json;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::collections::HashMap;
use straincompass_types::{GeneOrigin, OriginRecord};

/// Matches kept from the search; plenty to count plasmid vs chromosome.
const HITLIST: usize = 100;
/// Fewer full-length matches than this in the genus offers a wider search.
const MIN_MATCHES: usize = 5;
/// NCBI: poll any one request id no more often than once a minute.
const POLL_EVERY_SECS: i64 = 60;
/// A search NCBI has kept queued this long is given up on, so the user
/// gets an answer ("busy, try later") instead of an endless spinner.
const GIVE_UP_SECS: i64 = 45 * 60;
/// A match counts when it covers this much of the gene at this identity.
const MATCH_COV: f64 = 90.0;
const MATCH_PID: f64 = 90.0;
/// Complete plasmids offered for comparison; fetch_record takes ≤ 1 Mb.
const MAX_PLASMIDS: usize = 5;
const MAX_CHROMOSOMES: usize = 3;
const BACTERIA: &str = "all bacteria";

/// The stored state of one lookup, per gene sequence and genus.
#[derive(Serialize, Deserialize, Default)]
struct Stored {
    /// The search in flight, if any.
    #[serde(default)]
    rid: Option<String>,
    #[serde(default)]
    scope: String,
    #[serde(default)]
    scope_note: String,
    #[serde(default)]
    last_poll: i64,
    #[serde(default)]
    submitted: i64,
    #[serde(default)]
    result: Option<GeneOrigin>,
}

fn now() -> i64 {
    std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_secs() as i64)
        .unwrap_or(0)
}

/// One lookup at a time: two clicks must not submit two searches.
static LOCK: tokio::sync::Mutex<()> = tokio::sync::Mutex::const_new(());

#[derive(Deserialize)]
pub struct OriginQuery {
    pub gene_id: String,
    /// Search all bacteria instead of the project's genus. Explicit,
    /// because that search sits in NCBI's queue for many minutes.
    #[serde(default)]
    pub wide: bool,
}

/// GET /runs/{id}/panel_origin?gene_id : start or continue the lookup.
pub async fn panel_origin(
    State(state): State<SharedState>,
    Path(run_id): Path<i64>,
    Query(q): Query<OriginQuery>,
) -> ApiResult<Json<GeneOrigin>> {
    let (project_id, query_ids, status) = jobs::run_meta(&state, run_id)?;
    if status != "succeeded" {
        return Err(ApiError::BadRequest(
            "This run has not finished yet. Please wait for it to complete.".into(),
        ));
    }
    // a gene can hold several sequences; search with the one the strains
    // carry, not merely the first (for cadA that was an Enterococcus one)
    let variant = run_variant(&state, project_id, run_id, &query_ids, &q.gene_id);
    let seq = panel_gene_seq(&state, project_id, run_id, &variant)?;
    let searched_with = {
        let ref_genes: Vec<straincompass_types::WgaGene> =
            jobs::load_reference_json(&state, project_id, run_id)
                .ok()
                .and_then(|v| serde_json::from_value(v["genes"].clone()).ok())
                .unwrap_or_default();
        jobs::panel_record(&state, project_id, run_id, &variant)
            .map(|r| {
                straincompass_engine::element::panel_gene_source(&q.gene_id, &r.desc, &ref_genes)
            })
            .unwrap_or_default()
    };
    let genus = project_genus(&state, project_id);
    let (scope, scope_note) = match (&genus, q.wide) {
        (Some(g), false) => (
            g.clone(),
            format!("Searched {g} records (the project's genus)."),
        ),
        (None, false) => (
            BACTERIA.to_string(),
            "The project names no organism, so all bacteria were searched.".to_string(),
        ),
        (_, true) => (BACTERIA.to_string(), "Searched all bacteria.".to_string()),
    };
    let key = {
        let mut h = Sha256::new();
        h.update(&seq);
        h.update(b"|");
        h.update(scope.as_bytes());
        format!("{:x}", h.finalize())
    };
    let dir = state.data_dir.join("gene_origin");
    std::fs::create_dir_all(&dir)?;
    let path = dir.join(format!("{key}.json"));

    let _guard = LOCK.lock().await;
    let mut st: Stored = std::fs::read(&path)
        .ok()
        .and_then(|b| serde_json::from_slice(&b).ok())
        .unwrap_or_default();
    let save = |st: &Stored| {
        if let Ok(b) = serde_json::to_vec_pretty(st) {
            let _ = std::fs::write(&path, b);
        }
    };
    if let Some(mut r) = st.result.clone() {
        r.gene_id = q.gene_id.clone();
        r.searched_with = searched_with;
        return Ok(Json(r));
    }
    let c = client();
    let url = blast_url();
    let running = |st: &Stored, msg: &str| GeneOrigin {
        gene_id: q.gene_id.clone(),
        searched_with: searched_with.clone(),
        state: "running".into(),
        scope: st.scope.clone(),
        scope_note: st.scope_note.clone(),
        message: msg.into(),
        ..Default::default()
    };

    let Some(rid) = st.rid.clone() else {
        // nothing in flight: submit
        st.rid = Some(submit(&c, &url, &seq, &scope).await?);
        st.scope = scope.clone();
        st.scope_note = scope_note.clone();
        st.last_poll = now();
        st.submitted = now();
        save(&st);
        return Ok(Json(running(&st, "Search submitted to NCBI BLAST.")));
    };

    if st.submitted > 0 && now() - st.submitted > GIVE_UP_SECS {
        st.rid = None;
        save(&st);
        return Err(ApiError::Internal(
            "NCBI BLAST kept the search queued for over 45 minutes and is too busy right now. Please try again later; a new search will be started.".into(),
        ));
    }
    if now() - st.last_poll < POLL_EVERY_SECS {
        return Ok(Json(running(&st, "NCBI BLAST is still searching.")));
    }
    st.last_poll = now();
    let text = match poll(&c, &url, &rid).await {
        Ok(Some(t)) => t,
        Ok(None) => {
            save(&st);
            return Ok(Json(running(&st, "NCBI BLAST is still searching.")));
        }
        Err(e) => {
            // a lost search is resubmitted on the next request
            st.rid = None;
            save(&st);
            return Err(e);
        }
    };
    let qlen = seq.len() as u64;
    let matches = full_length_matches(&text, qlen);

    let summaries = summarize(&c, matches.keys().cloned().collect()).await?;
    let mut result = classify(&matches, &summaries);
    result.gene_id = q.gene_id.clone();
    result.searched_with = searched_with.clone();
    result.state = "done".into();
    result.scope = st.scope.clone();
    result.scope_note = st.scope_note.clone();
    if st.scope != BACTERIA && result.n_matches < MIN_MATCHES {
        // too little in the genus to say much; offer the slow search
        result.can_widen = true;
        result.scope_note.push_str(&format!(
            " Only {} full-length match{} there; searching all bacteria may say more.",
            result.n_matches,
            if result.n_matches == 1 { "" } else { "es" }
        ));
    }
    st.rid = None;
    st.result = Some(result.clone());
    save(&st);
    Ok(Json(result))
}

/// The panel record of `gene_id` that matched most often across the run:
/// full matches first, else any match; the gene's own id when nothing
/// matched or the run predates variants.
fn run_variant(
    state: &SharedState,
    project_id: i64,
    run_id: i64,
    query_ids: &[i64],
    gene_id: &str,
) -> String {
    use straincompass_types::Call;
    let mut present: HashMap<String, usize> = HashMap::new();
    let mut any: HashMap<String, usize> = HashMap::new();
    for qid in query_ids {
        let Ok(res) = jobs::load_query_result(state, project_id, run_id, *qid) else {
            continue;
        };
        let Some(r) = res
            .panel
            .unwrap_or_default()
            .into_iter()
            .find(|r| r.gene_id == gene_id && !r.variant.is_empty())
        else {
            continue;
        };
        if r.call == Call::Present {
            *present.entry(r.variant.clone()).or_default() += 1;
        }
        *any.entry(r.variant).or_default() += 1;
    }
    let pick = |m: HashMap<String, usize>| {
        m.into_iter()
            .max_by(|a, b| a.1.cmp(&b.1).then(b.0.cmp(&a.0)))
            .map(|(v, _)| v)
    };
    pick(present)
        .or_else(|| pick(any))
        .unwrap_or_else(|| gene_id.to_string())
}

/// The gene's sequence from the run's own panel copy.
fn panel_gene_seq(
    state: &SharedState,
    project_id: i64,
    run_id: i64,
    gene_id: &str,
) -> ApiResult<Vec<u8>> {
    let panel = state
        .run_dir(project_id, run_id)
        .join("panel")
        .join("panel.fa");
    let recs = straincompass_engine::fasta::parse_fasta(&panel).map_err(|_| {
        ApiError::NotFound("This run's gene panel file is no longer on the server.".into())
    })?;
    recs.into_iter()
        .find(|r| r.id == gene_id)
        .map(|r| r.seq)
        .filter(|s| s.len() >= 50)
        .ok_or_else(|| ApiError::NotFound(format!("{gene_id} is not in this run's panel.")))
}

/// The genus to search in: the project's organism, else the first word
/// of the reference FASTA description ("NC_003210.1 Listeria
/// monocytogenes EGD-e ..." -> "Listeria").
fn project_genus(state: &SharedState, project_id: i64) -> Option<String> {
    let organism: String = {
        let conn = state.db.lock().unwrap();
        conn.query_row(
            "SELECT organism FROM projects WHERE id = ?1",
            [project_id],
            |r| r.get::<_, Option<String>>(0),
        )
        .ok()
        .flatten()
        .unwrap_or_default()
    };
    if let Some(g) = genus_of(&organism) {
        return Some(g);
    }
    let (fasta, _) = crate::routes::uploads::reference_paths(state, project_id).ok()??;
    let mut first = String::new();
    use std::io::BufRead;
    std::io::BufReader::new(std::fs::File::open(fasta).ok()?)
        .read_line(&mut first)
        .ok()?;
    let desc = first
        .trim_start_matches('>')
        .split_once(char::is_whitespace)?
        .1;
    genus_of(desc)
}

/// A genus is a capitalised Latin word; anything else is not trusted as
/// an NCBI organism query.
fn genus_of(text: &str) -> Option<String> {
    let w = text.split_whitespace().next()?;
    let ok = w.len() >= 3
        && w.chars().next()?.is_ascii_uppercase()
        && w.chars().skip(1).all(|c| c.is_ascii_lowercase());
    ok.then(|| w.to_string())
}

async fn submit(c: &reqwest::Client, url: &str, seq: &[u8], scope: &str) -> ApiResult<String> {
    let entrez = if scope == BACTERIA {
        "Bacteria[Organism]".to_string()
    } else {
        format!("{scope}[Organism]")
    };
    let hitlist = HITLIST.to_string();
    let query = String::from_utf8_lossy(seq).into_owned();
    let resp = c
        .post(url)
        .form(&[
            ("CMD", "Put"),
            ("PROGRAM", "blastn"),
            ("MEGABLAST", "on"),
            ("DATABASE", "nt"),
            ("HITLIST_SIZE", hitlist.as_str()),
            ("ENTREZ_QUERY", entrez.as_str()),
            ("QUERY", query.as_str()),
        ])
        .timeout(std::time::Duration::from_secs(60))
        .send()
        .await
        .map_err(|e| ApiError::Internal(format!("NCBI BLAST could not be reached. ({e})")))?;
    let text = resp
        .text()
        .await
        .map_err(|e| ApiError::Internal(format!("NCBI BLAST answered oddly. ({e})")))?;
    text.lines()
        .find_map(|l| {
            l.trim()
                .strip_prefix("RID = ")
                .map(|r| r.trim().to_string())
        })
        .ok_or_else(|| ApiError::Internal("NCBI BLAST did not accept the search.".into()))
}

async fn poll(c: &reqwest::Client, url: &str, rid: &str) -> ApiResult<Option<String>> {
    let hitlist = HITLIST.to_string();
    let text = c
        .get(url)
        .query(&[
            ("CMD", "Get"),
            // FORMAT_TYPE=Tabular answers an empty page; the table comes
            // as text with the tabular alignment view
            ("FORMAT_TYPE", "Text"),
            ("ALIGNMENT_VIEW", "Tabular"),
            ("RID", rid),
            ("ALIGNMENTS", hitlist.as_str()),
            ("DESCRIPTIONS", hitlist.as_str()),
        ])
        .timeout(std::time::Duration::from_secs(60))
        .send()
        .await
        .map_err(|e| ApiError::Internal(format!("NCBI BLAST could not be reached. ({e})")))?
        .text()
        .await
        .map_err(|e| ApiError::Internal(format!("NCBI BLAST answered oddly. ({e})")))?;
    if text.contains("Status=UNKNOWN") {
        return Err(ApiError::Internal(
            "NCBI BLAST lost the search before it finished; it will be submitted again.".into(),
        ));
    }
    if text.contains("Status=WAITING") || text.contains("Status=QUEUED") {
        return Ok(None);
    }
    Ok(Some(text))
}

/// The accession in a BLAST subject id: "gi|123|gb|CP060527.1|" or a
/// bare "CP060527.1".
fn subject_accession(id: &str) -> String {
    let parts: Vec<&str> = id.split('|').filter(|p| !p.is_empty()).collect();
    parts
        .iter()
        .rev()
        .find(|p| p.contains('.') || p.chars().any(|c| c.is_ascii_alphabetic()) && p.len() > 4)
        .unwrap_or(&id)
        .to_string()
}

/// Full-length matches from a tabular result: accession -> (identity,
/// coverage of the gene), the best HSP per record.
fn full_length_matches(text: &str, qlen: u64) -> HashMap<String, (f64, f64)> {
    let mut out: HashMap<String, (f64, f64)> = HashMap::new();
    for line in text.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if line.starts_with('#') || f.len() < 8 {
            continue;
        }
        let (Ok(pid), Ok(qs), Ok(qe)) = (
            f[2].parse::<f64>(),
            f[6].parse::<u64>(),
            f[7].parse::<u64>(),
        ) else {
            continue;
        };
        let cov = 100.0 * (qe.max(qs) - qe.min(qs) + 1) as f64 / qlen.max(1) as f64;
        if cov < MATCH_COV || pid < MATCH_PID {
            continue;
        }
        let acc = subject_accession(f[1]);
        let e = out.entry(acc).or_insert((0.0, 0.0));
        if cov > e.1 || (cov == e.1 && pid > e.0) {
            *e = (pid, cov.min(100.0));
        }
    }
    out
}

/// Record summaries from NCBI: accession -> (title, length, genome type).
async fn summarize(
    c: &reqwest::Client,
    accessions: Vec<String>,
) -> ApiResult<HashMap<String, (String, u64, String)>> {
    let mut out = HashMap::new();
    for chunk in accessions.chunks(100) {
        let ids = chunk.join(",");
        let v: serde_json::Value = c
            .post("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi")
            .form(&[("db", "nuccore"), ("retmode", "json"), ("id", ids.as_str())])
            .timeout(std::time::Duration::from_secs(60))
            .send()
            .await
            .map_err(|e| ApiError::Internal(format!("NCBI could not be reached. ({e})")))?
            .json()
            .await
            .map_err(|e| ApiError::Internal(format!("NCBI answered oddly. ({e})")))?;
        if let Some(res) = v["result"].as_object() {
            for (uid, r) in res {
                if uid == "uids" {
                    continue;
                }
                let acc = r["accessionversion"]
                    .as_str()
                    .unwrap_or_default()
                    .to_string();
                out.insert(
                    acc,
                    (
                        r["title"].as_str().unwrap_or_default().to_string(),
                        r["slen"].as_u64().unwrap_or(0),
                        r["genome"].as_str().unwrap_or_default().to_string(),
                    ),
                );
            }
        }
    }
    Ok(out)
}

/// "plasmid", "chromosome" or "contig", from NCBI's genome type with the
/// title as a fallback.
fn record_kind(title: &str, genome: &str) -> &'static str {
    let t = title.to_lowercase();
    if genome == "plasmid" || t.contains(" plasmid") {
        "plasmid"
    } else if genome == "chromosome" || t.contains("chromosome") || t.contains("complete genome") {
        "chromosome"
    } else {
        "contig"
    }
}

fn classify(
    matches: &HashMap<String, (f64, f64)>,
    summaries: &HashMap<String, (String, u64, String)>,
) -> GeneOrigin {
    let mut g = GeneOrigin {
        n_matches: matches.len(),
        ..Default::default()
    };
    let mut plasmids = Vec::new();
    let mut chromosomes = Vec::new();
    for (acc, (pid, cov)) in matches {
        let (title, len, genome) = summaries.get(acc).cloned().unwrap_or_default();
        let kind = record_kind(&title, &genome);
        match kind {
            "plasmid" => g.n_plasmid += 1,
            "chromosome" => g.n_chromosome += 1,
            _ => g.n_contig += 1,
        }
        let rec = OriginRecord {
            accession: acc.clone(),
            title: title.clone(),
            length: len,
            kind: kind.into(),
            identity: *pid,
            coverage: *cov,
        };
        // a whole plasmid, not a gene record of one ("... plasmid pLM80
        // bcrB gene ..., complete CDS")
        let whole = title.to_lowercase().contains("complete sequence");
        if kind == "plasmid" && whole && (2_000..=1_000_000).contains(&len) {
            plasmids.push(rec);
        } else if kind == "chromosome" {
            chromosomes.push(rec);
        }
    }
    // closest copies first, the smaller plasmid on a tie: the comparison
    // reads best against a plasmid that is little more than the element
    plasmids.sort_by(|a, b| {
        b.identity
            .partial_cmp(&a.identity)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(a.length.cmp(&b.length))
    });
    chromosomes.sort_by(|a, b| {
        b.identity
            .partial_cmp(&a.identity)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    plasmids.truncate(MAX_PLASMIDS);
    chromosomes.truncate(MAX_CHROMOSOMES);
    g.plasmids = plasmids;
    g.chromosomes = chromosomes;
    g.message = if g.n_matches == 0 {
        "NCBI holds no record carrying this gene over its whole length.".into()
    } else {
        String::new()
    };
    g
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reads_genus_from_organism_or_header() {
        assert_eq!(
            genus_of("Listeria monocytogenes").as_deref(),
            Some("Listeria")
        );
        assert_eq!(
            genus_of("Bacillus cereus ATCC 14579 chromosome").as_deref(),
            Some("Bacillus")
        );
        assert_eq!(genus_of("").as_deref(), None);
        assert_eq!(genus_of("LM226_contig1").as_deref(), None);
        assert_eq!(genus_of("unknown bug").as_deref(), None);
    }

    #[test]
    fn takes_accessions_out_of_subject_ids() {
        assert_eq!(
            subject_accession("gi|1897675539|gb|CP060527.1|"),
            "CP060527.1"
        );
        assert_eq!(subject_accession("NZ_CP060527.1"), "NZ_CP060527.1");
        assert_eq!(subject_accession("ref|NC_003210.1|"), "NC_003210.1");
    }

    #[test]
    fn keeps_full_length_matches_best_hsp_per_record() {
        let text = "# BLASTN 2.16.0+\n\
            Query_1\tgb|CP060527.1|\t100.000\t924\t0\t0\t1\t924\t73097\t74020\t0.0\t1700\n\
            Query_1\tgb|CP060527.1|\t85.000\t100\t15\t0\t1\t100\t10\t110\t1e-20\t100\n\
            Query_1\tgb|CP000001.1|\t99.000\t500\t5\t0\t1\t500\t1\t500\t0.0\t900\n\
            Query_1\tgb|NC_000002.1|\t91.000\t900\t81\t0\t20\t919\t5\t904\t0.0\t1300\n";
        let m = full_length_matches(text, 924);
        assert_eq!(m.len(), 2, "{m:?}");
        assert_eq!(m["CP060527.1"], (100.0, 100.0));
        assert!((m["NC_000002.1"].1 - 97.40).abs() < 0.01);
    }

    #[test]
    fn parses_the_page_ncbi_actually_returns() {
        // FORMAT_TYPE=Text + ALIGNMENT_VIEW=Tabular, as received 2026-10-01
        let page = "<p><!--\nQBlastInfoBegin\n\tStatus=READY\nQBlastInfoEnd\n--><p>\n<PRE>\n\
            # blastn\n# Iteration: 0\n# Query: \n# RID: BW60W3C8016\n# Database: core_nt\n\
            # Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score\n\
            # 4 hits found\n\
            Query_3047295\tNG_076629.1\t100.000\t924\t0\t0\t1\t924\t101\t1024\t0.0\t1707\n\
            Query_3047295\tCP196590.1\t100.000\t924\t0\t0\t1\t924\t2047353\t2048276\t0.0\t1707\n\
            Query_3047295\tCP196591.1\t100.000\t924\t0\t0\t1\t924\t2323001\t2322078\t0.0\t1707\n\
            Query_3047295\tCP196592.1\t100.000\t924\t0\t0\t1\t924\t1265582\t1266505\t0.0\t1707\n\
            </PRE>\n";
        let m = full_length_matches(page, 924);
        assert_eq!(m.len(), 4, "{m:?}");
        assert_eq!(m["CP196591.1"], (100.0, 100.0));
    }

    #[test]
    fn sorts_records_into_plasmid_chromosome_and_contig() {
        assert_eq!(
            record_kind(
                "Listeria monocytogenes strain X plasmid pX, complete sequence",
                ""
            ),
            "plasmid"
        );
        assert_eq!(
            record_kind("Listeria monocytogenes EGD-e complete genome", ""),
            "chromosome"
        );
        assert_eq!(record_kind("whatever", "chromosome"), "chromosome");
        assert_eq!(
            record_kind(
                "Listeria monocytogenes strain Y NODE_1, whole genome shotgun sequence",
                ""
            ),
            "contig"
        );

        let mut matches = HashMap::new();
        matches.insert("P1".to_string(), (100.0, 100.0));
        matches.insert("P2".to_string(), (100.0, 100.0));
        matches.insert("C1".to_string(), (99.0, 100.0));
        matches.insert("W1".to_string(), (98.0, 100.0));
        matches.insert("G1".to_string(), (100.0, 100.0));
        let mut s = HashMap::new();
        s.insert(
            "P1".to_string(),
            (
                "x plasmid p1, complete sequence".to_string(),
                90_000,
                "plasmid".to_string(),
            ),
        );
        s.insert(
            "P2".to_string(),
            (
                "x plasmid p2, complete sequence".to_string(),
                50_000,
                "plasmid".to_string(),
            ),
        );
        s.insert(
            "G1".to_string(),
            (
                "x plasmid pLM80 bcrB gene, complete CDS".to_string(),
                518,
                "".to_string(),
            ),
        );
        s.insert(
            "C1".to_string(),
            (
                "x chromosome, complete genome".to_string(),
                3_000_000,
                "chromosome".to_string(),
            ),
        );
        s.insert(
            "W1".to_string(),
            (
                "x contig_5, whole genome shotgun sequence".to_string(),
                40_000,
                "".to_string(),
            ),
        );
        let g = classify(&matches, &s);
        // the gene record counts as a plasmid match but is not offered
        assert_eq!(
            (g.n_matches, g.n_plasmid, g.n_chromosome, g.n_contig),
            (5, 3, 1, 1)
        );
        let order: Vec<&str> = g.plasmids.iter().map(|p| p.accession.as_str()).collect();
        assert_eq!(order, vec!["P2", "P1"]);
        assert_eq!(g.chromosomes.len(), 1);
    }
}
