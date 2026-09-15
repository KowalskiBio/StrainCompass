//! Run lifecycle routes.

use crate::error::{ApiError, ApiResult};
use crate::models::{RunDto, RunLogDto, RunQueryDto};
use crate::state::SharedState;
use axum::extract::{Path, State};
use axum::Json;
use bactiment_types::{validate_params, RunParams};
use serde::Deserialize;
use std::sync::MutexGuard;

fn run_dto(
    conn: &MutexGuard<'_, rusqlite::Connection>,
    run_id: i64,
) -> ApiResult<Option<RunDto>> {
    let Ok(row) = conn.query_row(
        "SELECT id, project_id, status, step, error, created_at, started_at, finished_at, query_ids
         FROM runs WHERE id = ?1",
        [run_id],
        |r| {
            Ok((
                r.get::<_, i64>(0)?,
                r.get::<_, i64>(1)?,
                r.get::<_, String>(2)?,
                r.get::<_, Option<String>>(3)?,
                r.get::<_, Option<String>>(4)?,
                r.get::<_, String>(5)?,
                r.get::<_, Option<String>>(6)?,
                r.get::<_, Option<String>>(7)?,
                r.get::<_, String>(8)?,
            ))
        },
    ) else {
        return Ok(None);
    };
    let (id, project_id, status, step, error, created_at, started_at, finished_at, query_ids) = row;
    let has_panel: bool = conn
        .query_row(
            "SELECT EXISTS(SELECT 1 FROM files WHERE project_id = ?1 AND role = 'panel')",
            [project_id],
            |r| r.get(0),
        )
        .unwrap_or(false);
    let ids: Vec<i64> = serde_json::from_str(&query_ids).unwrap_or_default();
    let mut queries = Vec::new();
    for qid in ids {
        let name: Option<String> = conn
            .query_row(
                "SELECT display_name FROM files WHERE id = ?1",
                [qid],
                |r| r.get(0),
            )
            .ok();
        if let Some(name) = name {
            queries.push(RunQueryDto { file_id: qid, name });
        }
    }
    Ok(Some(RunDto {
        id,
        project_id,
        status,
        step,
        error,
        created_at,
        started_at,
        finished_at,
        queries,
        has_panel,
    }))
}

#[derive(Deserialize)]
pub struct StartRun {
    pub query_ids: Vec<i64>,
    #[serde(default)]
    pub params: Option<RunParams>,
}

/// POST /projects/{id}/runs
pub async fn start(
    State(state): State<SharedState>,
    Path(project_id): Path<i64>,
    Json(body): Json<StartRun>,
) -> ApiResult<Json<RunDto>> {
    crate::routes::uploads::ensure_project(&state, project_id).await?;
    if body.query_ids.is_empty() {
        return Err(ApiError::BadRequest(
            "Please choose at least one query genome to compare.".into(),
        ));
    }
    if body.query_ids.len() > 50 {
        return Err(ApiError::BadRequest(
            "Please compare at most 50 query genomes in one run.".into(),
        ));
    }
    let params = body.params.unwrap_or_default();
    if let Err(errors) = validate_params(&params) {
        return Err(ApiError::BadRequest(
            serde_json::to_string(&errors).unwrap_or_else(|_| "Invalid parameters.".into()),
        ));
    }
    // inputs must exist and belong to the project
    {
        let conn = state.db.lock().unwrap();
        for qid in &body.query_ids {
            let ok: bool = conn
                .query_row(
                    "SELECT EXISTS(SELECT 1 FROM files WHERE id = ?1 AND project_id = ?2 AND role = 'query')",
                    rusqlite::params![qid, project_id],
                    |r| r.get(0),
                )
                .unwrap_or(false);
            if !ok {
                return Err(ApiError::BadRequest(
                    "One of the selected query genomes does not exist anymore.".into(),
                ));
            }
        }
        // reference must be present
        let has_ref: bool = conn
            .query_row(
                "SELECT EXISTS(SELECT 1 FROM files WHERE project_id = ?1 AND role = 'reference_fasta')",
                [project_id],
                |r| r.get(0),
            )
            .unwrap_or(false);
        if !has_ref {
            return Err(ApiError::BadRequest(
                "Please add a reference genome (FASTA + GFF) before comparing.".into(),
            ));
        }
    }
    let run_id = {
        let conn = state.db.lock().unwrap();
        conn.execute(
            "INSERT INTO runs (project_id, status, params_json, query_ids) VALUES (?1, 'queued', ?2, ?3)",
            rusqlite::params![
                project_id,
                serde_json::to_string(&params).unwrap(),
                serde_json::to_string(&body.query_ids).unwrap(),
            ],
        )?;
        conn.last_insert_rowid()
    };
    crate::jobs::spawn_run(state.clone(), run_id);
    let conn = state.db.lock().unwrap();
    let dto = run_dto(&conn, run_id)?.ok_or_else(|| ApiError::Internal("run vanished".into()))?;
    Ok(Json(dto))
}

/// GET /runs/{id}
pub async fn detail(State(state): State<SharedState>, Path(run_id): Path<i64>) -> ApiResult<Json<RunLogDto>> {
    let (dto, logs) = {
        let conn = state.db.lock().unwrap();
        let dto = run_dto(&conn, run_id)?
            .ok_or_else(|| ApiError::NotFound("This run does not exist (anymore).".into()))?;
        let mut stmt = conn.prepare("SELECT line FROM run_logs WHERE run_id = ?1 ORDER BY seq")?;
        let logs = stmt
            .query_map([run_id], |r| r.get::<_, String>(0))?
            .filter_map(|l| l.ok())
            .collect::<Vec<_>>();
        (dto, logs)
    };
    Ok(Json(RunLogDto { run: dto, logs }))
}

/// GET /projects/{id}/runs
pub async fn list_for_project(
    State(state): State<SharedState>,
    Path(project_id): Path<i64>,
) -> ApiResult<Json<Vec<RunDto>>> {
    let conn = state.db.lock().unwrap();
    let mut stmt = conn
        .prepare("SELECT id FROM runs WHERE project_id = ?1 ORDER BY id DESC")?;
    let ids: Vec<i64> = stmt
        .query_map([project_id], |r| r.get(0))?
        .collect::<std::result::Result<_, _>>()?;
    drop(stmt);
    let mut out = Vec::new();
    for id in ids {
        if let Some(dto) = run_dto(&conn, id)? {
            out.push(dto);
        }
    }
    Ok(Json(out))
}

/// DELETE /runs/{id}
pub async fn delete(State(state): State<SharedState>, Path(run_id): Path<i64>) -> ApiResult<&'static str> {
    let project_id: Option<i64> = {
        let conn = state.db.lock().unwrap();
        conn.query_row("SELECT project_id FROM runs WHERE id = ?1", [run_id], |r| r.get(0))
            .ok()
    };
    let Some(project_id) = project_id else {
        return Err(ApiError::NotFound("This run does not exist (anymore).".into()));
    };
    {
        let conn = state.db.lock().unwrap();
        conn.execute("DELETE FROM runs WHERE id = ?1", [run_id])?;
        conn.execute("DELETE FROM run_logs WHERE run_id = ?1", [run_id])?;
    }
    let dir = state.run_dir(project_id, run_id);
    let _ = std::fs::remove_dir_all(dir);
    Ok("deleted")
}

/// GET /runs/{id}/params : the parameter set used, plus the defaults.
pub async fn params(State(state): State<SharedState>, Path(run_id): Path<i64>) -> ApiResult<Json<serde_json::Value>> {
    let stored: Option<String> = {
        let conn = state.db.lock().unwrap();
        conn.query_row("SELECT params_json FROM runs WHERE id = ?1", [run_id], |r| r.get(0))
            .ok()
    };
    let Some(stored) = stored else {
        return Err(ApiError::NotFound("This run does not exist (anymore).".into()));
    };
    let params: RunParams = serde_json::from_str(&stored)
        .map_err(|_| ApiError::Internal("The stored parameters could not be read.".into()))?;
    Ok(Json(serde_json::json!({
        "params": params,
        "defaults": RunParams::default(),
        "schema": bactiment_types::param_schema(),
        "presets": ["default", "strict", "loose"],
    })))
}
