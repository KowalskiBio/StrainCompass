//! Locating and invoking the external tools (MUMmer, BLAST+).

use crate::{friendly, EngineError, Result};
use std::path::{Path, PathBuf};
use std::process::Command;

#[derive(Debug, Clone)]
pub struct ToolPaths {
    pub nucmer: PathBuf,
    pub show_coords: PathBuf,
    pub show_snps: PathBuf,
    pub dnadiff: PathBuf,
    pub makeblastdb: PathBuf,
    pub blastn: PathBuf,
    /// The gene finder, used only to predict the genes inside gained
    /// regions. Optional on purpose: `discover` is all-or-nothing, and
    /// every deployment made before gained regions existed lacks this
    /// binary, so requiring it would fail every run on every one of them.
    pub prodigal: Option<PathBuf>,
}

fn find_tool(name: &str, extra_dirs: &[PathBuf]) -> std::result::Result<PathBuf, EngineError> {
    for d in extra_dirs {
        let p = d.join(name);
        if p.is_file() {
            return Ok(p);
        }
    }
    if let Ok(path) = which(name) {
        return Ok(path);
    }
    Err(EngineError::ToolMissing(format!(
        "The tool {name} was not found. Install MUMmer / BLAST+ or set STRAINCOMPASS_TOOLS_DIRS."
    )))
}

fn which(name: &str) -> std::result::Result<PathBuf, EngineError> {
    let out = Command::new("which")
        .arg(name)
        .output()
        .map_err(|e| EngineError::ToolMissing(format!("which {name}: {e}")))?;
    if out.status.success() {
        let s = String::from_utf8_lossy(&out.stdout).trim().to_string();
        if !s.is_empty() {
            return Ok(PathBuf::from(s));
        }
    }
    Err(EngineError::ToolMissing(format!(
        "The tool {name} was not found on this server."
    )))
}

impl ToolPaths {
    /// Discover tools from extra dirs (from STRAINCOMPASS_TOOLS_DIRS, colon
    /// separated) then from PATH.
    pub fn discover() -> Result<ToolPaths> {
        let extra: Vec<PathBuf> = std::env::var("STRAINCOMPASS_TOOLS_DIRS")
            .map(|v| {
                v.split(':')
                    .filter(|s| !s.is_empty())
                    .map(PathBuf::from)
                    .collect()
            })
            .unwrap_or_default();
        Ok(ToolPaths {
            nucmer: find_tool("nucmer", &extra)?,
            show_coords: find_tool("show-coords", &extra)?,
            show_snps: find_tool("show-snps", &extra)?,
            dnadiff: find_tool("dnadiff", &extra)?,
            makeblastdb: find_tool("makeblastdb", &extra)?,
            blastn: find_tool("blastn", &extra)?,
            prodigal: find_tool("prodigal", &extra).ok(),
        })
    }

    /// Run nucmer with the given parameters. If `cache_delta` exists and
    /// is non-empty, it is reused and nucmer is not run.
    #[allow(clippy::too_many_arguments)]
    pub fn run_nucmer(
        &self,
        ref_fasta: &Path,
        qry_fasta: &Path,
        out_dir: &Path,
        prefix: &str,
        minmatch: Option<i64>,
        breaklen: Option<i64>,
        cache_delta: Option<&Path>,
    ) -> Result<(PathBuf, bool)> {
        let delta = out_dir.join(format!("{prefix}.delta"));
        if let Some(c) = cache_delta {
            if c.is_file() && std::fs::metadata(c).map(|m| m.len() > 0).unwrap_or(false) {
                std::fs::copy(c, &delta)?;
                return Ok((delta, true));
            }
        }
        let mut cmd = Command::new(&self.nucmer);
        // no matcher flag: nucmer default (--mumreference) matches the R
        // pipeline; --mum would lose alignments in query-repetitive regions
        cmd.arg("-p")
            .arg(out_dir.join(prefix))
            .arg(ref_fasta)
            .arg(qry_fasta);
        if let Some(l) = minmatch {
            cmd.arg("-l").arg(l.to_string());
        }
        if let Some(b) = breaklen {
            cmd.arg("-b").arg(b.to_string());
        }
        let out = cmd
            .output()
            .map_err(|e| EngineError::ToolMissing(format!("nucmer: {e}")))?;
        if !out.status.success() || !delta.is_file() {
            let msg = String::from_utf8_lossy(&out.stderr).trim().to_string();
            return Err(friendly(format!(
                "The genome alignment (nucmer) did not finish correctly. {msg}"
            )));
        }
        if let Some(c) = cache_delta {
            let _ = std::fs::copy(&delta, c);
        }
        Ok((delta, false))
    }

    /// Run dnadiff for the overall report.
    pub fn run_dnadiff(
        &self,
        ref_fasta: &Path,
        qry_fasta: &Path,
        out_dir: &Path,
        prefix: &str,
    ) -> Result<PathBuf> {
        let report = out_dir.join(format!("{prefix}.report"));
        let out = Command::new(&self.dnadiff)
            .arg("-p")
            .arg(out_dir.join(prefix))
            .arg(ref_fasta)
            .arg(qry_fasta)
            .output()
            .map_err(|e| EngineError::ToolMissing(format!("dnadiff: {e}")))?;
        if !out.status.success() || !report.is_file() {
            let msg = String::from_utf8_lossy(&out.stderr).trim().to_string();
            return Err(friendly(format!(
                "The overall comparison report (dnadiff) could not be computed. {msg}"
            )));
        }
        // dnadiff always writes <prefix>.snps, a full SNP listing that runs
        // ~19 MB per bacterial genome pair. Nothing reads it: not the engine,
        // not the API, and it is not in the Files panel whitelist. It was 305
        // MB of a 440 MB project on disk. Drop it once the report exists;
        // failing to delete it must never fail the run.
        let _ = std::fs::remove_file(out_dir.join(format!("{prefix}.snps")));
        Ok(report)
    }
}
