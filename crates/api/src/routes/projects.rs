//! Project CRUD and usage.

use crate::error::{ApiError, ApiResult};
use crate::models::ProjectDto;
use crate::state::SharedState;
use axum::extract::{Path, State};
use axum::Json;
use serde::Deserialize;
use std::sync::MutexGuard;

pub fn all_projects(
    conn: &MutexGuard<'_, rusqlite::Connection>,
    state: &crate::state::AppState,
) -> ApiResult<Vec<ProjectDto>> {
    let mut dtos = Vec::new();
    let mut stmt = conn.prepare("SELECT id FROM projects ORDER BY created_at DESC, id DESC")?;
    let ids: Vec<i64> = stmt
        .query_map([], |r| r.get(0))?
        .collect::<std::result::Result<_, _>>()?;
    drop(stmt);
    for id in ids {
        if let Some(dto) = one(conn, state, id)? {
            dtos.push(dto);
        }
    }
    Ok(dtos)
}

#[derive(Deserialize)]
pub struct CreateProject {
    pub name: String,
    pub organism: Option<String>,
}

pub async fn create(
    State(state): State<SharedState>,
    Json(body): Json<CreateProject>,
) -> ApiResult<Json<ProjectDto>> {
    let name = body.name.trim().to_string();
    if name.is_empty() {
        return Err(ApiError::BadRequest(
            "Please give the project a name.".into(),
        ));
    }
    if name.len() > 200 {
        return Err(ApiError::BadRequest(
            "The project name is too long (more than 200 characters).".into(),
        ));
    }
    let organism = body
        .organism
        .as_deref()
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .unwrap_or("bacteria")
        .to_string();
    let conn = state.db.lock().unwrap();
    conn.execute(
        "INSERT INTO projects (name, organism) VALUES (?1, ?2)",
        (&name, &organism),
    )?;
    let id = conn.last_insert_rowid();
    drop(conn);
    std::fs::create_dir_all(state.project_dir(id))?;
    let dto = ProjectDto {
        id,
        name,
        organism,
        created_at: crate::db::now_rfc3339(),
        n_runs: 0,
        n_queries: 0,
        has_reference: false,
        has_panel: false,
        usage_bytes: 0,
    };
    Ok(Json(dto))
}

pub async fn list(State(state): State<SharedState>) -> ApiResult<Json<Vec<ProjectDto>>> {
    let conn = state.db.lock().unwrap();
    let dtos = all_projects(&conn, &state)?;
    drop(conn);
    Ok(Json(dtos))
}

pub fn one(
    conn: &MutexGuard<'_, rusqlite::Connection>,
    state: &crate::state::AppState,
    id: i64,
) -> ApiResult<Option<ProjectDto>> {
    let Ok(row) = conn.query_row(
        "SELECT id, name, organism, created_at FROM projects WHERE id = ?1",
        [id],
        |r| {
            Ok((
                r.get::<_, i64>(0)?,
                r.get::<_, String>(1)?,
                r.get::<_, String>(2)?,
                r.get::<_, String>(3)?,
            ))
        },
    ) else {
        return Ok(None);
    };
    let (id, name, organism, created_at) = row;
    let n_runs: i64 = conn.query_row(
        "SELECT COUNT(*) FROM runs WHERE project_id = ?1",
        [id],
        |r| r.get(0),
    )?;
    let n_queries: i64 = conn.query_row(
        "SELECT COUNT(*) FROM files WHERE project_id = ?1 AND role = 'query'",
        [id],
        |r| r.get(0),
    )?;
    let has_reference: i64 = conn.query_row(
        "SELECT EXISTS(SELECT 1 FROM files WHERE project_id = ?1 AND role = 'reference_fasta')",
        [id],
        |r| r.get(0),
    )?;
    let has_panel: i64 = conn.query_row(
        "SELECT EXISTS(SELECT 1 FROM files WHERE project_id = ?1 AND role = 'panel')",
        [id],
        |r| r.get(0),
    )?;
    let usage = crate::files::dir_size(&state.project_dir(id));
    Ok(Some(ProjectDto {
        id,
        name,
        organism,
        created_at,
        n_runs,
        n_queries,
        has_reference: has_reference != 0,
        has_panel: has_panel != 0,
        usage_bytes: usage,
    }))
}

pub async fn detail(
    State(state): State<SharedState>,
    Path(id): Path<i64>,
) -> ApiResult<Json<ProjectDto>> {
    let conn = state.db.lock().unwrap();
    let dto = one(&conn, &state, id)?
        .ok_or_else(|| ApiError::NotFound("This project does not exist (anymore).".into()))?;
    drop(conn);
    Ok(Json(dto))
}

#[derive(Deserialize)]
pub struct RenameProject {
    pub name: String,
}

pub async fn rename(
    State(state): State<SharedState>,
    Path(id): Path<i64>,
    Json(body): Json<RenameProject>,
) -> ApiResult<Json<ProjectDto>> {
    let name = body.name.trim().to_string();
    if name.is_empty() {
        return Err(ApiError::BadRequest(
            "Please give the project a name.".into(),
        ));
    }
    {
        let conn = state.db.lock().unwrap();
        let n = conn.execute(
            "UPDATE projects SET name = ?1 WHERE id = ?2",
            rusqlite::params![name, id],
        )?;
        if n == 0 {
            return Err(ApiError::NotFound(
                "This project does not exist (anymore).".into(),
            ));
        }
    }
    detail(State(state), Path(id)).await
}

pub async fn delete(
    State(state): State<SharedState>,
    Path(id): Path<i64>,
) -> ApiResult<&'static str> {
    let conn = state.db.lock().unwrap();
    let n = conn.execute("DELETE FROM projects WHERE id = ?1", [id])?;
    if n == 0 {
        return Err(ApiError::NotFound(
            "This project does not exist (anymore).".into(),
        ));
    }
    drop(conn);
    let dir = state.project_dir(id);
    if dir.exists() {
        std::fs::remove_dir_all(&dir).map_err(|e| {
            ApiError::Internal(format!(
                "The project files could not be removed from the server. ({e})"
            ))
        })?;
    }
    Ok("deleted")
}
