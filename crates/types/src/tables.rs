//! Row types for the result tables and the viewer endpoints.

use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
#[serde(rename_all = "UPPERCASE")]
pub enum Call {
    #[default]
    Present,
    Partial,
    Absent,
}

impl Call {
    pub fn as_str(&self) -> &'static str {
        match self {
            Call::Present => "PRESENT",
            Call::Partial => "PARTIAL",
            Call::Absent => "ABSENT",
        }
    }
}

/// One row of the genes coverage table (per query).
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct GeneCoverageRow {
    pub locus_tag: String,
    pub symbol: String,
    pub biotype: String,
    pub seqid: String,
    pub start: u64,
    pub end: u64,
    pub length: u64,
    pub cov_bp: u64,
    pub cov_pct: f64,
    pub call: Call,
    /// Best block identity over this gene (0-100), 0 when nothing aligned.
    pub best_identity: f64,
    pub mismatches: u64,
    pub indels: u64,
}

/// One row of the unaligned gaps table.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct GapRow {
    pub seqid: String,
    pub start: u64,
    pub end: u64,
    pub length: u64,
    pub n_genes: usize,
    /// Locus tags of genes that overlap the gap.
    pub genes: Vec<String>,
}

/// One row of the strict panel recheck (per query).
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct PanelRow {
    pub gene_id: String,
    pub qlen: u64,
    pub cov_pct: f64,
    pub identity: f64,
    pub best_evalue: String,
    pub call: Call,
}

/// Presence/absence across all queries of a run.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct MatrixRow {
    pub locus_tag: String,
    pub symbol: String,
    pub biotype: String,
    pub seqid: String,
    pub start: u64,
    pub end: u64,
    /// One entry per query, same order as the run's query list.
    pub calls: Vec<Call>,
    pub cov_pcts: Vec<f64>,
}

/// Alignment block on the reference, for the genome viewer.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WgaBlock {
    pub ref_seqid: String,
    pub ref_start: u64,
    pub ref_end: u64,
    pub qry_seqid: String,
    pub qry_start: u64,
    pub qry_end: u64,
    pub qry_rev: bool,
    /// 0-100
    pub identity: f64,
}

/// Reference gene for the genome viewer track.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WgaGene {
    pub locus_tag: String,
    pub symbol: String,
    pub biotype: String,
    pub seqid: String,
    pub start: u64,
    pub end: u64,
    pub strand: i8,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WgaQuery {
    pub query_id: i64,
    pub query_name: String,
    pub blocks: Vec<WgaBlock>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WgaData {
    /// Reference sequence lengths.
    pub reference: Vec<(String, u64)>,
    pub genes: Vec<WgaGene>,
    pub queries: Vec<WgaQuery>,
}

/// Hover preview + alignment rows for one gene (MSA viewer).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeneBlock {
    /// Reference coordinates of the aligned slice (1-based inclusive).
    pub ref_start: u64,
    pub ref_end: u64,
    /// Query coordinates of the aligned slice (on the forward query
    /// sequence; original orientation preserved through qry_rev).
    pub qry_start: u64,
    pub qry_end: u64,
    pub qry_rev: bool,
    pub identity: f64,
    pub ref_seq: String,
    pub qry_seq: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeneQueryAlignment {
    pub query_id: i64,
    pub query_name: String,
    pub call: Call,
    pub cov_pct: f64,
    pub best_identity: f64,
    pub mismatches: u64,
    pub indels: u64,
    /// Alignment blocks overlapping the gene, ordered by reference position.
    pub blocks: Vec<GeneBlock>,
    /// Unaligned reference stretches between blocks, (start, end) inclusive.
    pub unaligned: Vec<(u64, u64)>,
    /// Premature stop codons found in the query sequence of this gene.
    pub premature_stops: Vec<PrematureStop>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PrematureStop {
    /// Codon index (0-based from the gene start).
    pub codon_index: u64,
    /// Amino acid position (1-based) where the stop appears.
    pub aa_position: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeneDetail {
    pub locus_tag: String,
    pub symbol: String,
    pub biotype: String,
    pub seqid: String,
    pub start: u64,
    pub end: u64,
    pub strand: i8,
    pub length: u64,
    pub reference_seq: String,
    pub queries: Vec<GeneQueryAlignment>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum RunStatus {
    Queued,
    Running,
    Succeeded,
    Failed,
}

impl RunStatus {
    pub fn as_str(&self) -> &'static str {
        match self {
            RunStatus::Queued => "queued",
            RunStatus::Running => "running",
            RunStatus::Succeeded => "succeeded",
            RunStatus::Failed => "failed",
        }
    }
}

/// Table sorting / filtering / pagination query, shared by the table
/// endpoints and the export endpoints so exports match the visible view.
#[derive(Debug, Clone, Deserialize)]
pub struct TableQuery {
    #[serde(default)]
    pub query_id: Option<i64>,
    #[serde(default = "default_page")]
    pub page: u64,
    #[serde(default = "default_page_size")]
    pub page_size: u64,
    #[serde(default)]
    pub sort_by: Option<String>,
    #[serde(default)]
    pub sort_dir: Option<String>,
    #[serde(default)]
    pub search: Option<String>,
    /// Filter by call: present | partial | absent.
    #[serde(default)]
    pub call: Option<String>,
    /// For exports: comma separated column keys; empty = all.
    #[serde(default)]
    pub cols: Option<String>,
}

fn default_page() -> u64 {
    0
}
fn default_page_size() -> u64 {
    200
}

/// A page of table rows plus total count, so the UI can paginate.
#[derive(Debug, Clone, Serialize)]
pub struct Page<T> {
    pub rows: Vec<T>,
    pub total: u64,
}
