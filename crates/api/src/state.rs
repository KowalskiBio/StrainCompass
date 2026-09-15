//! Application state: configuration, database, shared paths.

use std::path::PathBuf;
use std::sync::{Arc, Mutex};

pub struct AppState {
    pub db: Mutex<rusqlite::Connection>,
    pub data_dir: PathBuf,
    /// Limits how many comparisons run at once across all runs.
    pub cpu_slots: std::sync::Arc<tokio::sync::Semaphore>,
}

pub type SharedState = Arc<AppState>;

impl AppState {
    pub fn project_dir(&self, project_id: i64) -> PathBuf {
        self.data_dir.join("projects").join(project_id.to_string())
    }
    pub fn run_dir(&self, project_id: i64, run_id: i64) -> PathBuf {
        self.project_dir(project_id)
            .join("runs")
            .join(run_id.to_string())
    }
    pub fn uploads_dir(&self, project_id: i64) -> PathBuf {
        self.project_dir(project_id).join("uploads")
    }
    pub fn cache_dir(&self, project_id: i64) -> PathBuf {
        self.project_dir(project_id).join("cache")
    }
}
