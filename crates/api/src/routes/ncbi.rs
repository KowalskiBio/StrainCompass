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
        .user_agent("straincompass/0.1")
        .build()
        .expect("reqwest client")
}

/// An assembly accession: "GCF_000196035.1" or "GCA_000196035.1" (GCF
/// RefSeq, GCA GenBank), nine digits, optionally a version.
fn valid_accession(acc: &str) -> bool {
    let acc = acc.trim();
    let Some(rest) = acc
        .strip_prefix("GCF_")
        .or_else(|| acc.strip_prefix("GCA_"))
    else {
        return false;
    };
    let (digits, version) = match rest.split_once('.') {
        Some((d, v)) => (d, Some(v)),
        None => (rest, None),
    };
    digits.len() == 9
        && digits.bytes().all(|b| b.is_ascii_digit())
        && version.is_none_or(|v| !v.is_empty() && v.bytes().all(|b| b.is_ascii_digit()))
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
    // the local reference library first: its genomes are the files NCBI
    // serves, with no download
    if let Some((title, r)) = crate::routes::library::find_reference(&state, &accession) {
        let read = |p: &std::path::Path| {
            std::fs::read(p).map_err(|e| {
                ApiError::Internal(format!("The reference library could not be read. ({e})"))
            })
        };
        let fna = gunzip(&read(&r.fna_gz)?)?;
        let gff = gunzip(&read(&r.gff_gz)?)?;
        let name = if r.assembly_name.is_empty() {
            r.accession.clone()
        } else {
            r.assembly_name.replace([' ', ','], "_")
        };
        crate::routes::uploads::delete_role(&state, project_id, "reference_fasta").await?;
        crate::routes::uploads::delete_role(&state, project_id, "reference_gff").await?;
        let f1 = crate::routes::uploads::store_upload_public(
            &state,
            project_id,
            "reference_fasta",
            &format!("{} {name} genome.fna", r.accession),
            fna,
        )
        .await?;
        let f2 = crate::routes::uploads::store_upload_public(
            &state,
            project_id,
            "reference_gff",
            &format!("{} {name} annotation.gff", r.accession),
            gff,
        )
        .await?;
        return Ok(Json(
            serde_json::json!({ "fasta": f1, "gff": f2, "source": title }),
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
    // the note tells the run where the sequence came from
    let record = format!(">{gene} NCBI {acc}:{start}-{stop}\n{}", wrap(&seq));
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
    // "L28104" / "CP038643" / WGS "DBJORO010000002" style: 1-6 letters
    // then digits
    let plain = |c: &str| {
        let letters = c.chars().take_while(|c| c.is_ascii_uppercase()).count();
        let rest = &c[letters..];
        (1..=6).contains(&letters) && rest.len() >= 5 && rest.chars().all(|c| c.is_ascii_digit())
    };
    // RefSeq "NG_076629" / "NZ_CP168866" / "NZ_JAARYM010000001": a
    // two-letter prefix, then digits or a plain accession
    if let Some((pfx, num)) = core.split_once('_') {
        return pfx.len() == 2
            && pfx.chars().all(|c| c.is_ascii_uppercase())
            && ((num.len() >= 5 && num.chars().all(|c| c.is_ascii_digit())) || plain(num));
    }
    plain(core)
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
        record: format!(">{name} NCBI {}\n{}", spec.accession, wrap(&seq)),
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

/// Largest record `fetch_record` accepts: bigger than any plasmid, far
/// smaller than a chromosome worth comparing this way.
pub const MAX_RECORD_BP: usize = 1_000_000;

/// Fetch one whole nucleotide record (e.g. a complete plasmid) by
/// accession: (title, uppercase sequence).
pub async fn fetch_record(accession: &str, api_key: Option<&str>) -> ApiResult<(String, Vec<u8>)> {
    let acc = accession.trim();
    if !looks_like_accession(acc) {
        return Err(ApiError::BadRequest(format!(
            "\u{201c}{acc}\u{201d} is not a GenBank accession. Use the accession of a complete record, e.g. NZ_CP168866.1."
        )));
    }
    let mut url = format!(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={acc}&rettype=fasta&retmode=text"
    );
    if let Some(k) = api_key {
        url.push_str(&format!("&api_key={k}"));
    }
    let bytes = download(&client(), &url).await?;
    let text = String::from_utf8_lossy(&bytes);
    let mut lines = text.lines();
    let title = match lines.next() {
        Some(h) if h.starts_with('>') => h[1..].trim().to_string(),
        _ => {
            return Err(ApiError::BadRequest(format!(
                "NCBI has no nucleotide record {acc}."
            )))
        }
    };
    let seq: Vec<u8> = lines
        .take_while(|l| !l.starts_with('>'))
        .flat_map(|l| l.trim().bytes())
        .map(|b| b.to_ascii_uppercase())
        .collect();
    if seq.len() > MAX_RECORD_BP {
        return Err(ApiError::BadRequest(format!(
            "{acc} is {} bp long; only records up to 1 Mb (plasmids, transposons) can be compared this way.",
            seq.len()
        )));
    }
    if seq.is_empty() {
        return Err(ApiError::BadRequest(format!(
            "NCBI returned no sequence for {acc}."
        )));
    }
    Ok((title, seq))
}

/// A record title short enough for a note: "Listeria monocytogenes
/// transposon Tn5422 ATPase (cadA), accessory protein ..." keeps its
/// first 60 characters, cut at a word.
fn short_title(t: &str) -> String {
    if t.chars().count() <= 60 {
        return t.to_string();
    }
    let cut: String = t.chars().take(60).collect();
    let cut = cut.rsplit_once(' ').map(|(a, _)| a).unwrap_or(&cut);
    format!("{}\u{2026}", cut.trim_end_matches([',', ';']))
}

/// A record of the local reference library: the replicon, its sequence
/// and its annotated genes.
type LocalRecord = (
    straincompass_engine::library::Replicon,
    Vec<u8>,
    Vec<straincompass_engine::gff::Gene>,
);

fn local_record(
    lib: &straincompass_engine::library::Library,
    accession: &str,
) -> Option<LocalRecord> {
    let rep = lib.replicon(accession).ok()??;
    let seq = lib.replicon_seq(&rep).ok()??;
    let genes = lib.replicon_genes(&rep).ok()?;
    Some((rep, seq, genes))
}

/// How far around a panel record's origin to look for its operon
/// partners, bp each side: cadA and cadC of Tn5422 sit side by side, and
/// an operon rarely spans more than a few kb.
const PARTNER_WINDOW: u64 = 8_000;
/// At most this many records are scanned for partners per panel build,
/// to keep the NCBI traffic of one build small.
const MAX_PARTNER_RECORDS: usize = 15;

/// Variant sets for a freshly built panel: the genes beside each panel
/// record's origin that the panel holds under another source (the cadA
/// beside cadC in Tn5422), added as further variants, and a note for
/// every variant whose record is not from the project's genus. Returns
/// the records to append and the notes for the user. NCBI unreachable
/// means no additions, never a failed build.
///
/// A record the local reference library holds is read there, genes and
/// title alike; only the others go to NCBI.
pub async fn panel_variant_sets(
    fasta: &str,
    genus: &str,
    api_key: Option<&str>,
    lib: Option<&straincompass_engine::library::Library>,
) -> (String, Vec<String>) {
    use straincompass_engine::panel::variant_gene;
    use straincompass_engine::panel_variants::{
        parse_cds_fasta, partner_variants, record_origin, variant_summary, Neighbourhood,
    };
    let Ok(recs) = straincompass_engine::fasta::parse_fasta_str(fasta) else {
        return (String::new(), Vec::new());
    };
    let base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils";
    let delay = std::time::Duration::from_millis(if api_key.is_some() { 110 } else { 380 });
    let key = api_key.map(|k| format!("&api_key={k}")).unwrap_or_default();
    let c = client();
    let get = |url: String| {
        let c = c.clone();
        async move {
            tokio::time::sleep(delay).await;
            c.get(&url).send().await.ok()?.text().await.ok()
        }
    };

    let origins: Vec<(String, straincompass_engine::panel_variants::Origin)> = recs
        .iter()
        .filter_map(|r| record_origin(&r.desc).map(|o| (variant_gene(&r.id).to_string(), o)))
        .collect();
    // the records the library holds: (replicon, sequence, genes)
    let mut local: std::collections::HashMap<String, LocalRecord> = Default::default();
    if let Some(lib) = lib {
        for (_, o) in origins.iter().take(MAX_PARTNER_RECORDS) {
            if local.contains_key(&o.accession) {
                continue;
            }
            if let Some(rec) = local_record(lib, &o.accession) {
                local.insert(o.accession.clone(), rec);
            }
        }
    }
    let mut accessions: Vec<&str> = origins
        .iter()
        .map(|(_, o)| o.accession.as_str())
        .filter(|a| !local.contains_key(*a))
        .collect();
    accessions.sort();
    accessions.dedup();

    // record titles: what each origin is (organism, element)
    let mut titles: std::collections::HashMap<String, String> = local
        .iter()
        .map(|(acc, r)| (acc.clone(), r.0.title()))
        .collect();
    if !accessions.is_empty() {
        let url = format!(
            "{base}/esummary.fcgi?db=nuccore&id={}&retmode=json{key}",
            accessions.join(",")
        );
        if let Some(v) = get(url)
            .await
            .and_then(|t| serde_json::from_str::<serde_json::Value>(&t).ok())
        {
            if let Some(uids) = v["result"]["uids"].as_array() {
                for u in uids.iter().filter_map(|u| u.as_str()) {
                    let r = &v["result"][u];
                    if let (Some(acc), Some(t)) =
                        (r["accessionversion"].as_str(), r["title"].as_str())
                    {
                        titles.insert(acc.to_string(), short_title(t));
                    }
                }
            }
        }
    }

    let mut hoods: Vec<Neighbourhood> = Vec::new();
    let mut seen: Vec<straincompass_engine::panel_variants::Origin> = Vec::new();
    for (gene, o) in origins.iter().take(MAX_PARTNER_RECORDS) {
        if seen.contains(o) {
            continue;
        }
        seen.push(o.clone());
        if let Some((rep, seq, genes)) = local.get(&o.accession) {
            let (lo, hi) = (
                o.lo.saturating_sub(PARTNER_WINDOW).max(1),
                o.hi + PARTNER_WINDOW,
            );
            let cds = genes
                .iter()
                .filter(|g| g.biotype == "protein_coding" && g.end >= lo && g.start <= hi)
                .filter_map(|g| {
                    let s = straincompass_engine::library::Replicon::slice(
                        seq,
                        g.start,
                        g.end,
                        g.strand < 0,
                    )?;
                    Some(straincompass_engine::panel_variants::Cds {
                        gene: straincompass_engine::panel_variants::base_gene_name(&g.symbol)
                            .to_string(),
                        // the accession the panel names, so the variant
                        // reads as a partner of its origin
                        accession: o.accession.clone(),
                        start: g.start,
                        end: g.end,
                        seq: s,
                    })
                })
                .collect();
            hoods.push(Neighbourhood {
                source_gene: gene.clone(),
                origin: o.clone(),
                title: rep.title(),
                source: "Library".into(),
                cds,
            });
            continue;
        }
        let url = format!(
            "{base}/efetch.fcgi?db=nuccore&id={}&seq_start={}&seq_stop={}&rettype=fasta_cds_na&retmode=text{key}",
            o.accession,
            o.lo.saturating_sub(PARTNER_WINDOW).max(1),
            o.hi + PARTNER_WINDOW
        );
        let Some(text) = get(url).await else {
            continue;
        };
        hoods.push(Neighbourhood {
            source_gene: gene.clone(),
            origin: o.clone(),
            title: titles
                .get(&o.accession)
                .cloned()
                .unwrap_or_else(|| o.accession.clone()),
            source: "NCBI".into(),
            cds: parse_cds_fasta(&text, &o.accession),
        });
    }

    let additions = partner_variants(&recs, &hoods);
    let mut records = String::new();
    let mut notes: Vec<String> = Vec::new();
    for a in &additions {
        records.push_str(&a.record);
        notes.push(a.note.clone());
    }
    // a sequence from another genus may be the gene in another form
    if !genus.is_empty() {
        let mut flagged: Vec<(String, String)> = Vec::new();
        for (gene, o) in &origins {
            let Some(t) = titles.get(&o.accession) else {
                continue;
            };
            let organism: String = t.split_whitespace().take(2).collect::<Vec<_>>().join(" ");
            let first = t.split_whitespace().next().unwrap_or("");
            if !first.eq_ignore_ascii_case(genus)
                && !flagged.contains(&(gene.clone(), o.accession.clone()))
            {
                flagged.push((gene.clone(), o.accession.clone()));
                notes.push(format!(
                    "{gene}: the sequence from {} comes from {organism}, not {genus}. A gene of the same name can differ a lot between genera; check the results say which variant matched.",
                    o.accession
                ));
            }
        }
    }
    let mut all = recs;
    if let Ok(more) = straincompass_engine::fasta::parse_fasta_str(&records) {
        all.extend(more);
    }
    notes.extend(variant_summary(&all));
    (records, notes)
}

#[cfg(test)]
mod accession_tests {
    use super::{looks_like_accession, valid_accession};

    #[test]
    fn accepts_assembly_accessions() {
        for a in [
            "GCF_000196035.1",
            "GCA_000196035.1",
            " GCF_054558205.1 ",
            "GCF_000196035",
        ] {
            assert!(valid_accession(a), "{a}");
        }
        for a in [
            "GCF_00019603.1",
            "GCF_000196035.",
            "GCX_000196035.1",
            "NC_003210.1",
            "",
            "GCF_000196035.1a",
        ] {
            assert!(!valid_accession(a), "{a}");
        }
    }

    #[test]
    fn accepts_genbank_refseq_and_wgs_accessions() {
        for a in [
            "L28104.1",
            "CP038643.1",
            "HF565366",
            "NG_076629.1",
            "NZ_CP168866.1",
            "NZ_DBJORO010000002.1",
            "DBJORO010000002.1",
            "NC_003210.1",
        ] {
            assert!(looks_like_accession(a), "{a}");
        }
    }

    #[test]
    fn rejects_gene_names_and_locus_tags() {
        for a in [
            "emrC",
            "lmo0444",
            "LM4B_02324",
            "LM6179_RS03640",
            "ACTATD_RS15010",
            "ACCESSION",
            "inlA",
            "lm4b_02329",
            "",
        ] {
            assert!(!looks_like_accession(a), "{a}");
        }
    }
}

#[cfg(test)]
mod variant_set_tests {
    use super::panel_variant_sets;

    /// Live NCBI: the catalog's cadA_Lm (Enterococcus) and cadC_Lm
    /// (Tn5422) as a panel; the Tn5422 cadA beside cadC must be added.
    /// Run with CAD_CATALOG=<the cad entries of AMR_CDS.fa>
    /// `cargo test -p straincompass-api -- --ignored live_`.
    #[tokio::test]
    #[ignore]
    async fn live_adds_the_tn5422_cada_beside_cadc() {
        let cat = std::fs::read(std::env::var("CAD_CATALOG").unwrap()).unwrap();
        let c = vec![straincompass_engine::catalog::Catalog::amrfinder(&cat)];
        let org = Some("Listeria monocytogenes");
        let mut fasta = String::new();
        for name in ["cadA", "cadC"] {
            for (i, g) in straincompass_engine::catalog::lookup_all(name, org, &c)
                .iter()
                .enumerate()
            {
                fasta.push_str(&format!(
                    ">{} {} {}: {} [{}]\n{}\n",
                    straincompass_engine::panel::variant_id(name, i + 1),
                    g.source,
                    g.symbol,
                    g.product,
                    g.origin,
                    String::from_utf8_lossy(&g.seq)
                ));
            }
        }
        let (extra, notes) = panel_variant_sets(&fasta, "Listeria", None, None).await;
        for n in &notes {
            eprintln!("NOTE {n}");
        }
        for l in extra.lines().filter(|l| l.starts_with('>')) {
            eprintln!("{l}");
        }
        assert!(extra.contains("L28104.1:158-2293"), "{extra}");
    }
}
