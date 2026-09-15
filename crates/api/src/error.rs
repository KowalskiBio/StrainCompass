//! API error type with plain language messages.

use axum::http::StatusCode;
use axum::response::{IntoResponse, Response};
use axum::Json;
use serde_json::json;

#[derive(Debug, thiserror::Error)]
pub enum ApiError {
    #[error("{0}")]
    BadRequest(String),
    #[error("{0}")]
    NotFound(String),
    #[error("{0}")]
    Internal(String),
    #[error(transparent)]
    Engine(#[from] bactiment_engine::EngineError),
    #[error(transparent)]
    Io(#[from] std::io::Error),
    #[error(transparent)]
    Db(#[from] rusqlite::Error),
}

impl IntoResponse for ApiError {
    fn into_response(self) -> Response {
        let (status, msg) = match &self {
            ApiError::BadRequest(m) => (StatusCode::BAD_REQUEST, m.clone()),
            ApiError::NotFound(m) => (StatusCode::NOT_FOUND, m.clone()),
            ApiError::Engine(e) => (StatusCode::BAD_REQUEST, e.to_string()),
            ApiError::Io(e) => (
                StatusCode::INTERNAL_SERVER_ERROR,
                format!("A file operation on the server failed. Please try again. ({e})"),
            ),
            ApiError::Db(e) => (
                StatusCode::INTERNAL_SERVER_ERROR,
                format!("The server database reported a problem. Please try again. ({e})"),
            ),
            ApiError::Internal(m) => (StatusCode::INTERNAL_SERVER_ERROR, m.clone()),
        };
        if status.is_server_error() {
            tracing::error!(%status, %msg, "request failed");
        }
        (status, Json(json!({ "error": msg }))).into_response()
    }
}

pub type ApiResult<T> = std::result::Result<T, ApiError>;
