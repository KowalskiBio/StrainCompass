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

fn gunzip(data: &[u8]) -> ApiResult<Vec<u8>> {
    let mut out = Vec::new();
    let mut decoder = flate2::read::GzDecoder::new(data);
    decoder.read_to_end(&mut out).map_err(|_| {
        ApiError::BadRequest("The downloaded genome file could not be unpacked.".into())
    })?;
    Ok(out)
}
