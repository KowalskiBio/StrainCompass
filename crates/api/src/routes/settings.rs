//! Settings: the NCBI API key (masked in responses).

use crate::error::ApiResult;
use crate::state::SharedState;
use axum::extract::State;
use axum::Json;
use serde::Deserialize;
use serde_json::json;

#[derive(Deserialize)]
pub struct KeyBody {
    pub key: String,
}

/// PUT /settings/ncbi_api_key
pub async fn put_key(
    State(state): State<SharedState>,
    Json(body): Json<KeyBody>,
) -> ApiResult<Json<serde_json::Value>> {
    let key = body.key.trim().to_string();
    if key.is_empty() {
        return Err(crate::error::ApiError::BadRequest(
            "The API key cannot be empty. To remove the key, use the remove button.".into(),
        ));
    }
    if key.len() > 200 || !key.chars().all(|c| c.is_ascii_alphanumeric()) {
        return Err(crate::error::ApiError::BadRequest(
            "This does not look like an NCBI API key. It is a short combination of letters and digits.".into(),
        ));
    }
    let conn = state.db.lock().unwrap();
    crate::db::set_setting(&conn, "ncbi_api_key", &key)?;
    drop(conn);
    let resp = json!({ "status": "ok", "masked": mask(&key) });
    Ok(Json(resp))
}

/// DELETE /settings/ncbi_api_key
pub async fn delete_key(State(state): State<SharedState>) -> ApiResult<Json<serde_json::Value>> {
    let conn = state.db.lock().unwrap();
    crate::db::delete_setting(&conn, "ncbi_api_key")?;
    drop(conn);
    Ok(Json(json!({ "status": "ok", "masked": null })))
}

/// GET /settings
pub async fn get(State(state): State<SharedState>) -> ApiResult<Json<serde_json::Value>> {
    let (masked, has) = {
        let conn = state.db.lock().unwrap();
        let key = crate::db::get_setting(&conn, "ncbi_api_key")?;
        let has = key.is_some();
        (key.map(|k| mask(&k)), has)
    };
    Ok(Json(json!({
        "ncbi_api_key": masked,
        "has_ncbi_api_key": has,
    })))
}

fn mask(key: &str) -> String {
    if key.len() <= 6 {
        "*".repeat(key.len())
    } else {
        format!("{}**{}", &key[..3], &key[key.len() - 3..])
    }
}
