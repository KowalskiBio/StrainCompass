//! bactiment-engine: the analysis pipeline as a library with no web or
//! database dependencies. Ports the logic of the original R script
//! (compare_genome_vs_reference.R) to Rust.

pub mod blast;
pub mod coverage;
pub mod delta;
pub mod fasta;
pub mod gaps;
pub mod gff;
pub mod msa;
pub mod panel;
pub mod pipeline;
pub mod tools;

/// Errors that carry a user-facing, plain language message.
#[derive(Debug, thiserror::Error)]
pub enum EngineError {
    #[error("{0}")]
    Friendly(String),
    #[error("A needed tool is not available: {0}")]
    ToolMissing(String),
    #[error("{0}")]
    Io(#[from] std::io::Error),
}

pub type Result<T> = std::result::Result<T, EngineError>;

pub fn friendly(msg: impl Into<String>) -> EngineError {
    EngineError::Friendly(msg.into())
}
