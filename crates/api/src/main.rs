//! bactiment API server: HTTP + job runner + static SPA serving.

mod db;
mod error;
mod files;
mod jobs;
mod models;

mod state;

use axum::extract::DefaultBodyLimit;
use axum::routing::{delete, get, post, put};
use axum::Router;
use state::{AppState, SharedState};
use std::sync::{Arc, Mutex};

fn routes() -> Router<SharedState> {
    Router::new()
        .route("/health", get(|| async { "{\"status\":\"ok\"}" }))
        .route(
            "/projects",
            post(routes::projects::create).get(routes::projects::list),
        )
        .route(
            "/projects/{id}",
            get(routes::projects::detail)
                .patch(routes::projects::rename)
                .delete(routes::projects::delete),
        )
        .route("/projects/{id}/files", get(routes::uploads::list_files))
        .route(
            "/projects/{id}/files/{file_id}",
            delete(routes::uploads::delete_file),
        )
        .route(
            "/projects/{id}/reference",
            post(routes::uploads::upload_reference),
        )
        .route(
            "/projects/{id}/reference/ncbi",
            post(routes::ncbi::fetch_reference),
        )
        .route(
            "/projects/{id}/queries",
            post(routes::uploads::upload_queries),
        )
        .route("/projects/{id}/panel", post(routes::uploads::upload_panel))
        .route(
            "/projects/{id}/panel/from_ids",
            post(routes::uploads::upload_panel_ids),
        )
        .route(
            "/projects/{id}/panel/from_text",
            post(routes::uploads::upload_panel_text),
        )
        .route(
            "/projects/{id}/runs",
            post(routes::runs::start).get(routes::runs::list_for_project),
        )
        .route("/projects/{id}/usage", get(routes::usage::usage))
        .route(
            "/runs/{id}",
            get(routes::runs::detail).delete(routes::runs::delete),
        )
        .route("/runs/{id}/params", get(routes::runs::params))
        .route(
            "/runs/{id}/genes_coverage",
            get(routes::results::genes_coverage),
        )
        .route(
            "/runs/{id}/unaligned_gaps",
            get(routes::results::unaligned_gaps),
        )
        .route(
            "/runs/{id}/panel_recheck",
            get(routes::results::panel_recheck),
        )
        .route("/runs/{id}/matrix", get(routes::results::matrix))
        .route("/runs/{id}/wga", get(routes::results::wga))
        .route("/runs/{id}/alignment", get(routes::results::alignment))
        .route("/runs/{id}/refseq", get(routes::results::refseq))
        .route("/runs/{id}/gene/{locus}", get(routes::results::gene_detail))
        .route(
            "/runs/{id}/gene/{locus}/export",
            get(routes::runfiles::export_gene_alignment),
        )
        .route("/runs/{id}/files", get(routes::runfiles::list_run_files))
        .route(
            "/runs/{id}/files/{name}",
            get(routes::runfiles::get_run_file),
        )
        .route(
            "/runs/{id}/export/{table}",
            get(routes::export::export_table),
        )
        .route(
            "/settings/ncbi_api_key",
            put(routes::settings::put_key).delete(routes::settings::delete_key),
        )
        .route("/settings", get(routes::settings::get))
        .route("/presets", get(routes::presets::list))
}

fn spawn_maintenance() {
    // clean up stale runlogs
    tokio::spawn(async {
        let mut interval = tokio::time::interval(std::time::Duration::from_secs(600));
        loop {
            interval.tick().await;
            // placeholder for future maintenance
        }
    });
}

#[tokio::main]
async fn main() {
    tracing_subscriber::fmt()
        .with_env_filter(
            tracing_subscriber::EnvFilter::try_from_default_env()
                .unwrap_or_else(|_| "bactiment_api=info,tower_http=info".into()),
        )
        .init();

    let data_dir = std::env::var("BACTIMENT_DATA_DIR").unwrap_or_else(|_| "./data".into());
    let bind = std::env::var("BACTIMENT_BIND").unwrap_or_else(|_| "127.0.0.1:8010".into());
    let data_dir = std::path::PathBuf::from(data_dir);
    std::fs::create_dir_all(&data_dir).expect("cannot create the data directory");

    let conn = rusqlite::Connection::open(data_dir.join("bactiment.db"))
        .expect("cannot open the database");
    db::init_db(&conn).expect("cannot initialize the database");

    let cores = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(2);
    // keep one core for the web tier
    let slots = cores.saturating_sub(1).max(1);

    let state: SharedState = Arc::new(AppState {
        db: Mutex::new(conn),
        data_dir: data_dir.clone(),
        cpu_slots: std::sync::Arc::new(tokio::sync::Semaphore::new(slots)),
    });

    // SPA static serving: explicit BACTIMENT_STATIC_DIR, then frontend/dist
    // next to the data dir, then ./frontend/dist
    let dist = std::env::var_os("BACTIMENT_STATIC_DIR")
        .map(std::path::PathBuf::from)
        .filter(|p| p.exists())
        .or_else(|| {
            data_dir
                .parent()
                .map(|p| p.join("frontend").join("dist"))
                .filter(|p| p.exists())
        })
        .or_else(|| {
            let cwd_dist = std::path::PathBuf::from("frontend/dist");
            cwd_dist.exists().then_some(cwd_dist)
        });
    let dist = dist.unwrap_or_else(|| data_dir.join("static"));

    let app = Router::new()
        .nest(
            "/api",
            routes()
                .with_state(state.clone())
                .layer(DefaultBodyLimit::max(600 * 1024 * 1024))
                // gzip the big JSON payloads (the alignment viewer ships
                // every query's variant events at once)
                .layer(tower_http::compression::CompressionLayer::new()),
        )
        .fallback_service(
            tower_http::services::ServeDir::new(&dist)
                .append_index_html_on_directories(true)
                .fallback(tower_http::services::ServeFile::new(
                    dist.join("index.html"),
                )),
        )
        // the SPA shell must always be revalidated so a new deploy is
        // picked up on reload; hashed assets cache fine on their own
        .layer(axum::middleware::from_fn(
            |req: axum::extract::Request, next: axum::middleware::Next| async move {
                let resp = next.run(req).await;
                let html = resp
                    .headers()
                    .get(axum::http::header::CONTENT_TYPE)
                    .and_then(|v| v.to_str().ok())
                    .is_some_and(|v| v.starts_with("text/html"));
                if html {
                    let (mut parts, body) = resp.into_parts();
                    parts.headers.insert(
                        axum::http::header::CACHE_CONTROL,
                        axum::http::HeaderValue::from_static("no-cache"),
                    );
                    return axum::response::Response::from_parts(parts, body);
                }
                resp
            },
        ));

    let listener = tokio::net::TcpListener::bind(&bind)
        .await
        .expect("cannot bind");
    tracing::info!(
        "bactiment listening on http://{bind} (data: {}, {} parallel slots)",
        data_dir.display(),
        slots
    );
    spawn_maintenance();
    axum::serve(listener, app)
        .with_graceful_shutdown(async {
            let _ = tokio::signal::ctrl_c().await;
        })
        .await
        .expect("server error");
}

mod routes {
    pub mod export;
    pub mod ncbi;
    pub mod projects;
    pub mod results;
    pub mod runfiles;
    pub mod runs;
    pub mod settings;
    pub mod uploads;

    pub mod usage {
        use crate::error::ApiResult;
        use crate::state::SharedState;
        use axum::extract::{Path, State};
        use axum::Json;
        use serde_json::json;

        pub async fn usage(
            State(state): State<SharedState>,
            Path(project_id): Path<i64>,
        ) -> ApiResult<Json<serde_json::Value>> {
            let bytes = crate::files::dir_size(&state.project_dir(project_id));
            Ok(Json(json!({
                "bytes": bytes,
                "human": crate::files::human_size(bytes),
            })))
        }
    }

    pub mod presets {
        use axum::Json;
        use serde_json::json;

        pub async fn list() -> Json<serde_json::Value> {
            Json(json!({
                "schema": bactiment_types::param_schema(),
                "presets": ["default", "strict", "loose"],
            }))
        }
    }
}
