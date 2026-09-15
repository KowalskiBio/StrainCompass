//! API response models.

use serde::Serialize;

#[derive(Debug, Clone, Serialize)]
pub struct ProjectDto {
    pub id: i64,
    pub name: String,
    pub organism: String,
    pub created_at: String,
    pub n_runs: i64,
    pub n_queries: i64,
    pub has_reference: bool,
    pub has_panel: bool,
    pub usage_bytes: u64,
}

#[derive(Debug, Clone, Serialize)]
pub struct FileDto {
    pub id: i64,
    pub role: String,
    pub display_name: String,
    pub size: u64,
    pub created_at: String,
}

#[derive(Debug, Clone, Serialize)]
pub struct RunDto {
    pub id: i64,
    pub project_id: i64,
    pub status: String,
    pub step: Option<String>,
    pub error: Option<String>,
    pub created_at: String,
    pub started_at: Option<String>,
    pub finished_at: Option<String>,
    pub queries: Vec<RunQueryDto>,
    pub has_panel: bool,
}

#[derive(Debug, Clone, Serialize)]
pub struct RunQueryDto {
    pub file_id: i64,
    pub name: String,
}

#[derive(Debug, Clone, Serialize)]
pub struct RunLogDto {
    pub run: RunDto,
    pub logs: Vec<String>,
}
