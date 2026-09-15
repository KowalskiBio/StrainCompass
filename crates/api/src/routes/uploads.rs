//! Upload routes: reference (files or NCBI accession), queries, panel.

use crate::error::{ApiError, ApiResult};
use crate::models::FileDto;
use crate::state::SharedState;
use axum::extract::{Multipart, Path, State};
use axum::Json;
use bactiment_engine::fasta;
use std::io::Write;
use uuid::Uuid;

const MAX_UPLOAD_BYTES: usize = 512 * 1024 * 1024;

fn safe_filename(name: &str) -> String {
    let cleaned: String = name
        .chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || matches!(c, '.' | '-' | '_') {
                c
            } else {
                '_'
            }
        })
        .collect();
    if cleaned.is_empty() {
        "file".into()
    } else {
        cleaned
    }
}

/// Validate + store one upload. Returns the DB row.
pub async fn store_upload_public(
    state: &SharedState,
    project_id: i64,
    role: &str,
    display_name: &str,
    bytes: Vec<u8>,
) -> ApiResult<FileDto> {
    store_upload(state, project_id, role, display_name, bytes).await
}

async fn store_upload(
    state: &SharedState,
    project_id: i64,
    role: &str,
    display_name: &str,
    bytes: Vec<u8>,
) -> ApiResult<FileDto> {
    if bytes.len() > MAX_UPLOAD_BYTES {
        return Err(ApiError::BadRequest(
            "This file is too large. The limit is 512 MB.".into(),
        ));
    }
    if bytes.is_empty() {
        return Err(ApiError::BadRequest(
            "The uploaded file is empty. Please check the file and try again.".into(),
        ));
    }
    // Content validation with friendly errors.
    let text = String::from_utf8_lossy(&bytes);
    match role {
        "reference_fasta" | "query" | "panel" => {
            fasta::parse_fasta_str(&text).map_err(|e| {
                ApiError::BadRequest(format!(
                    "\u{201c}{display_name}\u{201d}: {}",
                    e
                ))
            })?;
        }
        "reference_gff" => {
            bactiment_engine::gff::parse_gff_str(&text).map_err(|e| {
                ApiError::BadRequest(format!("\u{201c}{display_name}\u{201d}: {}", e))
            })?;
        }
        _ => {}
    }
    let uploads = state.uploads_dir(project_id);
    std::fs::create_dir_all(&uploads)?;
    let stored_name = format!("{}_{}", Uuid::new_v4().simple(), safe_filename(display_name));
    let path = uploads.join(&stored_name);
    let mut f = std::fs::File::create(&path)?;
    f.write_all(&bytes)?;
    drop(f);
    let created_at = crate::db::now_rfc3339();
    let id = {
        let conn = state.db.lock().unwrap();
        conn.execute(
            "INSERT INTO files (project_id, role, display_name, stored_name, size, created_at)
             VALUES (?1, ?2, ?3, ?4, ?5, ?6)",
            rusqlite::params![project_id, role, display_name, stored_name, bytes.len() as i64, created_at],
        )?;
        conn.last_insert_rowid()
    };
    Ok(FileDto {
        id,
        role: role.into(),
        display_name: display_name.to_string(),
        size: bytes.len() as u64,
        created_at,
    })
}

/// POST /projects/{id}/reference : multipart fields `fasta`, `gff`.
pub async fn upload_reference(
    State(state): State<SharedState>,
    Path(project_id): Path<i64>,
    mut multipart: Multipart,
) -> ApiResult<Json<Vec<FileDto>>> {
    ensure_project(&state, project_id).await?;
    let mut fasta_bytes: Option<(String, Vec<u8>)> = None;
    let mut gff_bytes: Option<(String, Vec<u8>)> = None;
    while let Some(field) = multipart
        .next_field()
        .await
        .map_err(|_| ApiError::BadRequest("The upload could not be read. Please try again.".into()))?
    {
        let name = field.name().unwrap_or("").to_string();
        let filename = field.file_name().unwrap_or("file").to_string();
        let data = field
            .bytes()
            .await
            .map_err(|_| ApiError::BadRequest("The upload was interrupted. Please try again.".into()))?;
        match name.as_str() {
            "fasta" => fasta_bytes = Some((filename, data.to_vec())),
            "gff" => gff_bytes = Some((filename, data.to_vec())),
            _ => {}
        }
    }
    let (fname, fdata) = fasta_bytes.ok_or_else(|| {
        ApiError::BadRequest(
            "Please provide both the reference genome (FASTA) and its annotation (GFF).".into(),
        )
    })?;
    let (gname, gdata) = gff_bytes.ok_or_else(|| {
        ApiError::BadRequest(
            "Please provide both the reference genome (FASTA) and its annotation (GFF).".into(),
        )
    })?;
    let mut created = Vec::new();
    {
        // replace any previous reference
        delete_role(&state, project_id, "reference_fasta").await?;
        delete_role(&state, project_id, "reference_gff").await?;
    }
    created.push(store_upload(&state, project_id, "reference_fasta", &fname, fdata).await?);
    created.push(store_upload(&state, project_id, "reference_gff", &gname, gdata).await?);
    Ok(Json(created))
}

/// POST /projects/{id}/queries : one or many FASTA files.
pub async fn upload_queries(
    State(state): State<SharedState>,
    Path(project_id): Path<i64>,
    mut multipart: Multipart,
) -> ApiResult<Json<Vec<FileDto>>> {
    ensure_project(&state, project_id).await?;
    let mut created = Vec::new();
    while let Some(field) = multipart
        .next_field()
        .await
        .map_err(|_| ApiError::BadRequest("The upload could not be read. Please try again.".into()))?
    {
        let filename = field.file_name().unwrap_or("").to_string();
        let data = field
            .bytes()
            .await
            .map_err(|_| ApiError::BadRequest("The upload was interrupted. Please try again.".into()))?;
        if filename.is_empty() {
            continue;
        }
        created.push(
            store_upload(&state, project_id, "query", &filename, data.to_vec()).await?,
        );
    }
    if created.is_empty() {
        return Err(ApiError::BadRequest(
            "No query files were received. Please choose at least one FASTA file.".into(),
        ));
    }
    Ok(Json(created))
}

/// POST /projects/{id}/panel : one FASTA file.
pub async fn upload_panel(
    State(state): State<SharedState>,
    Path(project_id): Path<i64>,
    mut multipart: Multipart,
) -> ApiResult<Json<FileDto>> {
    ensure_project(&state, project_id).await?;
    let mut got: Option<(String, Vec<u8>)> = None;
    while let Some(field) = multipart
        .next_field()
        .await
        .map_err(|_| ApiError::BadRequest("The upload could not be read. Please try again.".into()))?
    {
        let filename = field.file_name().unwrap_or("").to_string();
        let data = field
            .bytes()
            .await
            .map_err(|_| ApiError::BadRequest("The upload was interrupted. Please try again.".into()))?;
        if filename.is_empty() {
            continue;
        }
        got = Some((filename, data.to_vec()));
    }
    let (name, data) = got.ok_or_else(|| {
        ApiError::BadRequest("No gene panel file was received. Please choose a FASTA file.".into())
    })?;
    delete_role(&state, project_id, "panel").await?;
    let dto = store_upload(&state, project_id, "panel", &name, data).await?;
    Ok(Json(dto))
}

/// GET /projects/{id}/files
pub async fn list_files(
    State(state): State<SharedState>,
    Path(project_id): Path<i64>,
) -> ApiResult<Json<Vec<FileDto>>> {
    ensure_project(&state, project_id).await?;
    let conn = state.db.lock().unwrap();
    let mut stmt = conn.prepare(
        "SELECT id, role, display_name, size, created_at FROM files
         WHERE project_id = ?1 ORDER BY role, created_at",
    )?;
    let rows = stmt
        .query_map([project_id], |r| {
            Ok(FileDto {
                id: r.get(0)?,
                role: r.get(1)?,
                display_name: r.get(2)?,
                size: r.get::<_, i64>(3)? as u64,
                created_at: r.get(4)?,
            })
        })?
        .collect::<std::result::Result<Vec<_>, _>>()?;
    Ok(Json(rows))
}

/// DELETE /projects/{id}/files/{file_id}
pub async fn delete_file(
    State(state): State<SharedState>,
    Path((project_id, file_id)): Path<(i64, i64)>,
) -> ApiResult<&'static str> {
    ensure_project(&state, project_id).await?;
    let stored_name: Option<String> = {
        let conn = state.db.lock().unwrap();
        conn.query_row(
            "SELECT stored_name FROM files WHERE id = ?1 AND project_id = ?2",
            [file_id, project_id],
            |r| r.get(0),
        )
        .ok()
    };
    let Some(stored_name) = stored_name else {
        return Err(ApiError::NotFound("This file does not exist (anymore).".into()));
    };
    {
        let conn = state.db.lock().unwrap();
        conn.execute("DELETE FROM files WHERE id = ?1", [file_id])?;
    }
    let path = state.uploads_dir(project_id).join(stored_name);
    let _ = std::fs::remove_file(path);
    Ok("deleted")
}

pub async fn delete_role(state: &SharedState, project_id: i64, role: &str) -> ApiResult<()> {
    let names: Vec<String> = {
        let conn = state.db.lock().unwrap();
        let mut stmt =
            conn.prepare("SELECT stored_name FROM files WHERE project_id = ?1 AND role = ?2")?;
        let rows = stmt
            .query_map(rusqlite::params![project_id, role], |r| r.get(0))?
            .collect::<std::result::Result<Vec<_>, _>>()?;
        rows
    };
    if !names.is_empty() {
        let conn = state.db.lock().unwrap();
        conn.execute(
            "DELETE FROM files WHERE project_id = ?1 AND role = ?2",
            rusqlite::params![project_id, role],
        )?;
        drop(conn);
        for n in names {
            let _ = std::fs::remove_file(state.uploads_dir(project_id).join(n));
        }
    }
    Ok(())
}

pub async fn ensure_project(state: &SharedState, project_id: i64) -> ApiResult<()> {
    let exists = {
        let conn = state.db.lock().unwrap();
        conn.query_row("SELECT 1 FROM projects WHERE id = ?1", [project_id], |_| Ok(()))
            .is_ok()
    };
    if !exists {
        return Err(ApiError::NotFound("This project does not exist (anymore).".into()));
    }
    Ok(())
}

/// Resolve the current reference (fasta, gff) of a project.
pub fn reference_paths(
    state: &SharedState,
    project_id: i64,
) -> ApiResult<Option<(std::path::PathBuf, std::path::PathBuf)>> {
    let conn = state.db.lock().unwrap();
    let mut stmt = conn.prepare(
        "SELECT role, stored_name FROM files WHERE project_id = ?1 AND role IN ('reference_fasta','reference_gff')",
    )?;
    let rows = stmt
        .query_map([project_id], |r| {
            Ok((r.get::<_, String>(0)?, r.get::<_, String>(1)?))
        })?
        .collect::<std::result::Result<Vec<_>, _>>()?;
    drop(stmt);
    let mut fasta = None;
    let mut gff = None;
    for (role, stored) in rows {
        match role.as_str() {
            "reference_fasta" => fasta = Some(state.uploads_dir(project_id).join(stored)),
            "reference_gff" => gff = Some(state.uploads_dir(project_id).join(stored)),
            _ => {}
        }
    }
    Ok(match (fasta, gff) {
        (Some(f), Some(g)) => Some((f, g)),
        _ => None,
    })
}
