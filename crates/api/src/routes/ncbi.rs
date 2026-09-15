//! Fetch a reference genome (FASTA + GFF) from NCBI by accession.

use crate::error::{ApiError, ApiResult};
use crate::routes::uploads::ensure_project;
use crate::state::SharedState;
use axum::extract::{Path, State};
use axum::Json;
use serde::Deserialize;
use std::io::Read;
use std::time::Duration;

fn client() -> reqwest::Client {
    reqwest::Client::builder()
        .timeout(Duration::from_secs(120))
        .connect_timeout(Duration::from_secs(20))
        .user_agent("bactiment/0.1")
        .build()
        .expect("reqwest client")
}

fn valid_accession(acc: &str) -> bool {
    let re = |p: &str| {
        let mut it = p.split('_');
        match (it.next(), it.next(), it.next(), it.next()) {
            (Some(pfx), Some(v), Some(n), None) => {
                (pfx == "GCF" || pfx == "GCA")
                    && v.len() == 3
                    && v.chars().all(|c| c.is_ascii_digit())
                    && !n.is_empty()
                    && n.chars().all(|c| c.is_ascii_digit())
            }
            _ => false,
        }
    };
    let acc = acc.trim();
    re(acc)
}

#[derive(Deserialize)]
pub struct AccessionBody {
    pub accession: String,
}

/// POST /projects/{id}/reference/ncbi
/// Fetches the genome + annotation for the accession and stores them as
/// the project reference. Returns the created file entries.
pub async fn fetch_reference(
    State(state): State<SharedState>,
    Path(project_id): Path<i64>,
    Json(body): Json<AccessionBody>,
) -> ApiResult<Json<serde_json::Value>> {
    ensure_project(&state, project_id).await?;
    let accession = body.accession.trim().to_string();
    if !valid_accession(&accession) {
        return Err(ApiError::BadRequest(
            "This does not look like an NCBI assembly accession. It should look like GCF_000196035.1.".into(),
        ));
    }
    let api_key = {
        let conn = state.db.lock().unwrap();
        crate::db::get_setting(&conn, "ncbi_api_key")?
    };
    let c = client();

    // 1. esearch for the assembly uid
    let mut esearch = format!(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=assembly&term={}%5BAssembly%20Accession%5D&retmode=json",
        accession
    );
    if let Some(k) = &api_key {
        esearch.push_str(&format!("&api_key={k}"));
    }
    let resp = c.get(&esearch).send().await.map_err(|e| {
        ApiError::Internal(format!(
            "NCBI could not be reached. Please try again in a moment. ({e})"
        ))
    })?;
    if !resp.status().is_success() {
        return Err(ApiError::BadRequest(format!(
            "NCBI did not answer the search for {accession}. Please try again in a moment."
        )));
    }
    let json: serde_json::Value = resp.json().await.map_err(|_| {
        ApiError::BadRequest("NCBI returned an unexpected answer. Please try again.".into())
    })?;
    let uid = json["esearchresult"]["idlist"][0]
        .as_str()
        .ok_or_else(|| {
            ApiError::BadRequest(format!(
                "No assembly named {accession} was found on NCBI. Please check the accession."
            ))
        })?
        .to_string();

    // 2. esummary for the FTP path
    let mut esummary = format!(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=assembly&id={uid}&retmode=json"
    );
    if let Some(k) = &api_key {
        esummary.push_str(&format!("&api_key={k}"));
    }
    let resp = c.get(&esummary).send().await.map_err(|e| {
        ApiError::Internal(format!(
            "NCBI could not be reached. Please try again in a moment. ({e})"
        ))
    })?;
    let json: serde_json::Value = resp.json().await.map_err(|_| {
        ApiError::BadRequest("NCBI returned an unexpected answer. Please try again.".into())
    })?;
    let result = &json["result"][&uid];
    let ftp = result["ftppath_refseq"]
        .as_str()
        .filter(|s| !s.is_empty())
        .or_else(|| result["ftppath_genbank"].as_str().filter(|s| !s.is_empty()))
        .ok_or_else(|| {
            ApiError::BadRequest(format!(
                "No genome files are available for {accession} on NCBI."
            ))
        })?
        .to_string();
    let basename = ftp.rsplit('/').next().unwrap_or("").to_string();
    let display_name = result["assemblyname"]
        .as_str()
        .map(|s| s.replace([' ', ','], "_"))
        .unwrap_or_else(|| accession.clone());

    // 3. download fna.gz + gff.gz
    let fna_url = format!("{ftp}/{basename}_genomic.fna.gz");
    let gff_url = format!("{ftp}/{basename}_genomic.gff.gz");
    let fna_gz = download(&c, &fna_url).await?;
    let gff_gz = download(&c, &gff_url).await?;
    let fna = gunzip(&fna_gz)?;
    let gff = gunzip(&gff_gz)?;

    // 4. store as project reference
    let fasta_name = format!("{accession} {display_name} genome.fna");
    let gff_name = format!("{accession} {display_name} annotation.gff");
    crate::routes::uploads::delete_role(&state, project_id, "reference_fasta").await?;
    crate::routes::uploads::delete_role(&state, project_id, "reference_gff").await?;
    let f1 = crate::routes::uploads::store_upload_public(
        &state,
        project_id,
        "reference_fasta",
        &fasta_name,
        fna,
    )
    .await?;
    let f2 = crate::routes::uploads::store_upload_public(
        &state,
        project_id,
        "reference_gff",
        &gff_name,
        gff,
    )
    .await?;
    Ok(Json(serde_json::json!({ "fasta": f1, "gff": f2 })))
}

async fn download(c: &reqwest::Client, url: &str) -> ApiResult<Vec<u8>> {
    let resp = c.get(url).send().await.map_err(|e| {
        ApiError::Internal(format!(
            "The download from NCBI failed. Please try again. ({e})"
        ))
    })?;
    if !resp.status().is_success() {
        return Err(ApiError::BadRequest(
            "The genome files could not be downloaded from NCBI. Please try again in a moment."
                .into(),
        ));
    }
    let bytes = resp.bytes().await.map_err(|e| {
        ApiError::Internal(format!("The download from NCBI was interrupted. ({e})"))
    })?;
    if bytes.len() > 2 * MAX_DL {
        return Err(ApiError::BadRequest(
            "This genome is unexpectedly large (more than 1 GB).".into(),
        ));
    }
    Ok(bytes.to_vec())
}

const MAX_DL: usize = 512 * 1024 * 1024;

/// A gene sequence fetched from the NCBI Gene database.
pub struct NcbiGene {
    /// FASTA record (">name\nseq..."), header uses the user's spelling.
    pub record: String,
    /// Human-readable provenance, e.g. "NC_003212.1 (Listeria innocua)".
    pub source: String,
}

/// Look up one gene by name (symbol or locus tag) within a genus and
/// fetch its sequence from the annotated genome it sits on.
///
/// Returns None when NCBI has no such gene for this organism; ambiguous
/// hits prefer the project's own species.
pub async fn fetch_gene(
    c: &reqwest::Client,
    api_key: Option<&str>,
    organism: &str,
    genus: &str,
    gene: &str,
) -> Option<NcbiGene> {
    let base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils";
    // NCBI allows 3 requests/s without an API key, 10 with one.
    let delay = if api_key.is_some() { 110 } else { 380 };
    async fn eget(c: &reqwest::Client, url: String, delay_ms: u64) -> Option<reqwest::Response> {
        tokio::time::sleep(std::time::Duration::from_millis(delay_ms)).await;
        c.get(&url).send().await.ok()
    }
    let q = |s: &str| -> String { s.replace(' ', "+") };

    // 1. find gene ids, by [Gene Name] then by [All Fields]
    let mut uids: Vec<String> = Vec::new();
    for field in ["Gene%20Name", "All%20Fields"] {
        let term = format!(
            "{}%5B{}%5D%20AND%20{}%5BOrganism%5D",
            q(gene),
            field,
            q(genus)
        );
        let mut url = format!("{base}/esearch.fcgi?db=gene&term={term}&retmode=json&retmax=5");
        if let Some(k) = api_key {
            url.push_str(&format!("&api_key={k}"));
        }
        let resp = eget(c, url, delay).await?;
        let json: serde_json::Value = resp.json().await.ok()?;
        let list: Vec<String> = json["esearchresult"]["idlist"]
            .as_array()
            .map(|a| {
                a.iter()
                    .filter_map(|v| v.as_str().map(str::to_string))
                    .collect()
            })
            .unwrap_or_default();
        if !list.is_empty() {
            uids = list;
            break;
        }
    }
    if uids.is_empty() {
        return None;
    }

    // 2. summaries: order the hits (project species first, then any)
    let mut url = format!(
        "{base}/esummary.fcgi?db=gene&id={}&retmode=json",
        uids.join(",")
    );
    if let Some(k) = api_key {
        url.push_str(&format!("&api_key={k}"));
    }
    let resp = eget(c, url, delay).await?;
    let json: serde_json::Value = resp.json().await.ok()?;
    let result = json["result"].clone();
    let mut entries: Vec<(String, serde_json::Value)> = uids
        .iter()
        .map(|u| (u.clone(), result[u].clone()))
        .collect();
    entries.sort_by_key(|(_, r)| {
        let sci = r["organism"]["scientificname"].as_str().unwrap_or("");
        let has_gi = !r["genomicinfo"]
            .as_array()
            .cloned()
            .unwrap_or_default()
            .is_empty();
        (
            if sci == organism { 0 } else { 1 },
            if has_gi { 0 } else { 1 },
        )
    });

    // coordinates: esummary first, then the gene XML (some records,
    // e.g. WGS-derived ones, only carry coordinates in the XML)
    let mut coords: Option<(String, i64, i64, &'static str)> = None;
    for (_, r) in &entries {
        if let Some(gi) = r["genomicinfo"].as_array().and_then(|a| a.first()) {
            let acc = gi["chraccver"].as_str().unwrap_or("").to_string();
            if let (Some(a), Some(b)) = (gi["chrstart"].as_i64(), gi["chrstop"].as_i64()) {
                if !acc.is_empty() {
                    let strand = if a > b { "2" } else { "1" };
                    coords = Some((acc, a.min(b) + 1, a.max(b) + 1, strand));
                    break;
                }
            }
        }
    }
    if coords.is_none() {
        for (uid, _) in &entries {
            let mut url = format!("{base}/efetch.fcgi?db=gene&id={uid}&retmode=xml");
            if let Some(k) = api_key {
                url.push_str(&format!("&api_key={k}"));
            }
            let resp = eget(c, url, delay).await?;
            let xml = resp.text().await.ok()?;
            if let Some(xc) = parse_gene_xml(&xml) {
                coords = Some(xc);
                break;
            }
        }
    }
    let (acc, start, stop, strand) = coords?;
    let sci = entries
        .first()
        .and_then(|(_, r)| r["organism"]["scientificname"].as_str())
        .unwrap_or(genus)
        .to_string();

    // 3. fetch the sequence region
    let mut url = format!(
        "{base}/efetch.fcgi?db=nuccore&id={acc}&seq_start={start}&seq_stop={stop}&strand={strand}&rettype=fasta&retmode=text"
    );
    if let Some(k) = api_key {
        url.push_str(&format!("&api_key={k}"));
    }
    let resp = eget(c, url, delay).await?;
    let text = resp.text().await.ok()?;
    let mut lines = text.lines();
    let _header = lines.next()?;
    let seq: String = lines.collect::<Vec<_>>().join("");
    if seq.len() < 30 {
        return None;
    }
    let record = format!(">{gene}\n{}", wrap(&seq));
    Some(NcbiGene {
        record,
        source: format!("{acc} ({sci})"),
    })
}

fn wrap(seq: &str) -> String {
    let mut out = String::new();
    let bytes = seq.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        let end = (i + 60).min(bytes.len());
        out.push_str(&seq[i..end]);
        out.push('\n');
        i = end;
    }
    out
}

/// A gene pinned to a specific GenBank record, e.g. "qacH (HF565366.1)"
/// or "emrC (CP038643.1:1496-1882 rev)". The coordinates are needed when
/// the record does not name the gene in its feature table.
pub struct AccessionSpec {
    pub accession: String,
    pub name: Option<String>,
    pub coords: Option<(i64, i64, bool)>,
}

fn looks_like_accession(s: &str) -> bool {
    let core = s.split('.').next().unwrap_or(s);
    // "NG_076629" / "NZ_JAARYM010000001" style (prefix + digits)
    if let Some((pfx, num)) = core.split_once('_') {
        return !pfx.is_empty()
            && pfx.len() <= 2
            && pfx.chars().all(|c| c.is_ascii_uppercase())
            && num.len() >= 5
            && num.chars().all(|c| c.is_ascii_digit());
    }
    // "L28104" / "HF565366" / "CP038643" style (letters then digits)
    let letters: String = core
        .chars()
        .take_while(|c| c.is_ascii_uppercase())
        .collect();
    let rest = &core[letters.len()..];
    !letters.is_empty()
        && letters.len() <= 2
        && rest.len() >= 5
        && rest.chars().all(|c| c.is_ascii_digit())
}

/// Recognize "name (ACC)", "name (ACC:start-end [rev])" or a bare "ACC"
/// in one line of a gene list.
pub fn parse_accession_spec(line: &str) -> Option<AccessionSpec> {
    let line = line.trim();
    // parenthesized spec: "emrC (CP038643.1:1496-1882 rev)"
    if let Some(open) = line.find('(') {
        let close = line.rfind(')')?;
        let inner = &line[open + 1..close];
        let name = line[..open].trim();
        let (acc_part, coord_part) = match inner.split_once(':') {
            Some((a, c)) => (a.trim(), Some(c.trim())),
            None => (inner.trim(), None),
        };
        if looks_like_accession(acc_part) {
            let coords = coord_part.and_then(parse_coord_part);
            return Some(AccessionSpec {
                accession: acc_part.to_string(),
                name: (!name.is_empty()).then(|| name.to_string()),
                coords,
            });
        }
    }
    // bare accession (optionally with coordinates)
    let mut rest = line;
    let mut coords = None;
    if let Some((acc, c)) = line.split_once(':') {
        rest = acc;
        coords = parse_coord_part(c.trim());
    }
    let tok = rest.split_whitespace().next()?;
    if looks_like_accession(tok) {
        return Some(AccessionSpec {
            accession: tok.to_string(),
            name: None,
            coords,
        });
    }
    None
}

fn parse_coord_part(c: &str) -> Option<(i64, i64, bool)> {
    let rev = c.ends_with("rev") || c.ends_with("complement");
    let c = c
        .trim_end_matches("rev")
        .trim_end_matches("complement")
        .trim();
    let (a, b) = c.split_once("..").or_else(|| c.split_once('-'))?;
    let start: i64 = a.trim().parse().ok()?;
    let end: i64 = b.trim().parse().ok()?;
    Some((start, end, rev))
}

/// Fetch the sequence of a gene pinned to a GenBank record. When no
/// coordinates are given, the gene is located via the record's feature
/// table (falling back to the whole record when it is not named).
pub async fn fetch_gene_by_accession(
    c: &reqwest::Client,
    api_key: Option<&str>,
    spec: &AccessionSpec,
) -> Option<NcbiGene> {
    let base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils";
    let delay = if api_key.is_some() { 110 } else { 380 };
    async fn eget(c: &reqwest::Client, url: String, delay_ms: u64) -> Option<String> {
        tokio::time::sleep(std::time::Duration::from_millis(delay_ms)).await;
        c.get(&url).send().await.ok()?.text().await.ok()
    }
    let url = |params: String| {
        let mut u = format!("{base}/efetch.fcgi?db=nuccore&{params}");
        if let Some(k) = api_key {
            u.push_str(&format!("&api_key={k}"));
        }
        u
    };

    let (start, stop, minus) = match spec.coords {
        Some((s, e, rev)) => (Some(s), Some(e), rev),
        None => {
            if let Some(name) = &spec.name {
                let mut gb = eget(
                    c,
                    url(format!("id={}&rettype=gb&retmode=text", spec.accession)),
                    delay,
                )
                .await?;
                // RefSeq CON records carry no features; the GenBank
                // original (same id without the NZ_ prefix) does
                if !gb.contains("FEATURES") && spec.accession.starts_with("NZ_") {
                    gb = eget(
                        c,
                        url(format!(
                            "id={}&rettype=gb&retmode=text",
                            spec.accession.trim_start_matches("NZ_")
                        )),
                        delay,
                    )
                    .await?;
                }
                // not named in this record: take the whole record
                find_gene_location(&gb, name).unwrap_or((None, None, false))
            } else {
                (None, None, false)
            }
        }
    };

    let mut params = format!("id={}&rettype=fasta&retmode=text", spec.accession);
    if let (Some(s), Some(e)) = (start, stop) {
        params.push_str(&format!(
            "&seq_start={s}&seq_stop={e}&strand={}",
            if minus { "2" } else { "1" }
        ));
    }
    let text = eget(c, url(params), delay).await?;
    let mut lines = text.lines();
    let _header = lines.next()?;
    let seq: String = lines.collect::<Vec<_>>().join("");
    if seq.len() < 30 {
        return None;
    }
    // a whole record only makes sense as a panel gene when it is small
    // (a transposon or cassette); whole chromosomes are rejected
    if start.is_none() && seq.len() > 25_000 {
        return None;
    }
    let name = spec.name.clone().unwrap_or_else(|| spec.accession.clone());
    let source = match (start, stop) {
        (Some(s), Some(e)) => format!(
            "{}:{}-{}{}",
            spec.accession,
            s,
            e,
            if minus { " rev" } else { "" }
        ),
        _ => format!("{} (whole record, {} bp)", spec.accession, seq.len()),
    };
    Some(NcbiGene {
        record: format!(">{name}\n{}", wrap(&seq)),
        source,
    })
}

/// Find the location of the gene named `name` in a GenBank feature table.
fn find_gene_location(gb: &str, name: &str) -> Option<(Option<i64>, Option<i64>, bool)> {
    let lines: Vec<&str> = gb.lines().collect();
    let wanted = format!("/gene=\"{}\"", name);
    let lower = name.to_lowercase();
    let mut i = 0;
    while i < lines.len() {
        let l = lines[i];
        // a feature key line: 5 spaces, a key, spaces, then the location
        if l.starts_with("     ") && !l.starts_with("          ") {
            let body = &l[5..];
            let key = body.split_whitespace().next().unwrap_or("");
            if key == "gene" {
                let loc = body[key.len()..].trim().to_string();
                // scan the qualifiers of this feature
                let mut j = i + 1;
                let mut named = false;
                while j < lines.len()
                    && (lines[j].starts_with("          /") || lines[j].starts_with("          "))
                {
                    if lines[j].contains(&wanted)
                        || lines[j]
                            .to_lowercase()
                            .contains(&format!("/gene=\"{}\"", lower))
                    {
                        named = true;
                    }
                    j += 1;
                }
                if named {
                    let minus = loc.contains("complement");
                    let mut nums = [0i64; 2];
                    let mut k = 0;
                    let bytes = loc.as_bytes();
                    let mut x = 0;
                    while x < bytes.len() && k < 2 {
                        if bytes[x].is_ascii_digit() {
                            let s = x;
                            while x < bytes.len() && bytes[x].is_ascii_digit() {
                                x += 1;
                            }
                            nums[k] = loc[s..x].parse().unwrap_or(0);
                            k += 1;
                        } else {
                            x += 1;
                        }
                    }
                    if k == 2 && nums[0] > 0 && nums[1] > 0 {
                        return Some((
                            Some(nums[0].min(nums[1])),
                            Some(nums[0].max(nums[1])),
                            minus,
                        ));
                    }
                    return None;
                }
                i = j;
                continue;
            }
        }
        i += 1;
    }
    None
}

/// Extract the genomic coordinates from a gene-db efetch XML record.
/// Coordinates are 0-based inclusive, like the esummary ones.
fn parse_gene_xml(xml: &str) -> Option<(String, i64, i64, &'static str)> {
    // nucleotide accession (skip protein accessions like WP_...)
    let mut acc: Option<String> = None;
    for key in ["Gene-source_src-str1>", "Gene-commentary_accession>"] {
        let mut from = 0;
        while let Some(p) = xml[from..].find(key) {
            let s = from + p + key.len();
            let e = xml[s..].find('<').map(|x| s + x)?;
            let val = &xml[s..e];
            if !val.starts_with("WP_") && val.contains('_') {
                acc = Some(val.to_string());
                break;
            }
            from = e;
        }
        if acc.is_some() {
            break;
        }
    }
    let acc = acc?;
    let f = xml.find("<Seq-interval_from>")? + "<Seq-interval_from>".len();
    let fe = xml[f..].find('<').map(|x| f + x)?;
    let from: i64 = xml[f..fe].trim().parse().ok()?;
    let t = xml.find("<Seq-interval_to>")? + "<Seq-interval_to>".len();
    let te = xml[t..].find('<').map(|x| t + x)?;
    let to: i64 = xml[t..te].trim().parse().ok()?;
    let window_end = (te + 400).min(xml.len());
    let strand = if xml[te..window_end].contains("minus") {
        "2"
    } else {
        "1"
    };
    Some((acc, from.min(to) + 1, from.max(to) + 1, strand))
}

fn gunzip(data: &[u8]) -> ApiResult<Vec<u8>> {
    let mut out = Vec::new();
    let mut decoder = flate2::read::GzDecoder::new(data);
    decoder.read_to_end(&mut out).map_err(|_| {
        ApiError::BadRequest("The downloaded genome file could not be unpacked.".into())
    })?;
    Ok(out)
}
