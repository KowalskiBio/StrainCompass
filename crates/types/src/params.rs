//! Engine parameter model. Every threshold the original R script hardcoded
//! is a named parameter with a default, a range and a UI control.
//!
//! Parameters are layered: `align` parameters change the alignment itself
//! (they force a real nucmer re-run), `postprocess` parameters only affect
//! how a finished alignment is interpreted (changing them is cheap).

use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ParamLayer {
    Align,
    Postprocess,
}

/// Full parameter set for one analysis run.
#[derive(Debug, Clone, Serialize, Deserialize, garde::Validate)]
#[garde(context(()))]
pub struct RunParams {
    /// Minimum unaligned region reported (bp).
    #[garde(range(min = 0, max = 100_000))]
    pub min_gap: u64,
    /// % of a gene aligned to call it present.
    #[garde(range(min = 0.0, max = 100.0))]
    pub present_cov: f64,
    /// > partial_cov and < present_cov => PARTIAL, below => ABSENT.
    #[garde(range(min = 0.0, max = 100.0))]
    pub partial_cov: f64,
    /// % query coverage for a panel gene to be present (strict BLAST recheck).
    #[garde(range(min = 0.0, max = 100.0))]
    pub blast_cov: f64,
    /// % identity for a panel gene to be present.
    #[garde(range(min = 0.0, max = 100.0))]
    pub blast_pid: f64,
    /// BLAST e-value cutoff.
    #[garde(range(min = 1e-50, max = 10.0))]
    pub blast_evalue: f64,
    /// nucmer -l (minimal match length). None = tool default.
    #[garde(skip)]
    pub nucmer_minmatch: Option<i64>,
    /// nucmer -b (breakpoint distance). None = tool default.
    #[garde(skip)]
    pub nucmer_breaklen: Option<i64>,
    /// Run the overall alignment report (dnadiff) or skip it.
    #[garde(skip)]
    pub dnadiff: bool,
}

impl Default for RunParams {
    fn default() -> Self {
        Self {
            min_gap: 200,
            present_cov: 95.0,
            partial_cov: 1.0,
            blast_cov: 90.0,
            blast_pid: 90.0,
            blast_evalue: 1e-10,
            nucmer_minmatch: None,
            nucmer_breaklen: None,
            dnadiff: true,
        }
    }
}

/// The named presets offered in the UI.
pub fn preset(name: &str) -> Option<RunParams> {
    match name {
        "default" => Some(RunParams::default()),
        "strict" => Some(RunParams {
            present_cov: 99.0,
            min_gap: 100,
            blast_cov: 95.0,
            blast_pid: 95.0,
            blast_evalue: 1e-20,
            ..RunParams::default()
        }),
        "loose" => Some(RunParams {
            present_cov: 85.0,
            min_gap: 500,
            blast_cov: 80.0,
            blast_pid: 85.0,
            blast_evalue: 1e-5,
            ..RunParams::default()
        }),
        _ => None,
    }
}

/// Names and metadata of every parameter, used by the UI to build the
/// parameter panel. Ranges are enforced on the backend.
pub fn param_schema() -> Vec<ParamSpec> {
    vec![
        ParamSpec {
            name: "min_gap".into(),
            label: "Minimum unaligned region".into(),
            help: "Unaligned stretches of the reference shorter than this (in base pairs) are not reported.".into(),
            layer: ParamLayer::Postprocess,
            kind: ParamKind::Int {
                default: 200,
                min: 0,
                max: 100_000,
            },
        },
        ParamSpec {
            name: "present_cov".into(),
            label: "Present above coverage %".into(),
            help: "A gene counts as present when at least this % of its bases are aligned.".into(),
            layer: ParamLayer::Postprocess,
            kind: ParamKind::Float {
                default: 95.0,
                min: 0.0,
                max: 100.0,
            },
        },
        ParamSpec {
            name: "partial_cov".into(),
            label: "Partial above coverage %".into(),
            help: "A gene is PARTIAL when its coverage is above this % but below the present threshold. Below it, the gene is reported as absent.".into(),
            layer: ParamLayer::Postprocess,
            kind: ParamKind::Float {
                default: 1.0,
                min: 0.0,
                max: 100.0,
            },
        },
        ParamSpec {
            name: "blast_cov".into(),
            label: "Panel gene coverage %".into(),
            help: "A panel gene counts as present in the strict recheck when at least this % of its bases are found.".into(),
            layer: ParamLayer::Postprocess,
            kind: ParamKind::Float {
                default: 90.0,
                min: 0.0,
                max: 100.0,
            },
        },
        ParamSpec {
            name: "blast_pid".into(),
            label: "Panel gene identity %".into(),
            help: "A panel gene counts as present in the strict recheck when its best match is at least this % identical.".into(),
            layer: ParamLayer::Postprocess,
            kind: ParamKind::Float {
                default: 90.0,
                min: 0.0,
                max: 100.0,
            },
        },
        ParamSpec {
            name: "blast_evalue".into(),
            label: "Match significance cutoff".into(),
            help: "How strict the gene panel recheck search is. Smaller values mean fewer chance matches.".into(),
            layer: ParamLayer::Postprocess,
            kind: ParamKind::Float {
                default: 1e-10,
                min: 1e-50,
                max: 10.0,
            },
        },
        ParamSpec {
            name: "nucmer_minmatch".into(),
            label: "Minimal match length".into(),
            help: "Advanced: the exact-match seed length used by the aligner. Leave empty for the tool default. Changing this triggers a full re-alignment.".into(),
            layer: ParamLayer::Align,
            kind: ParamKind::OptionalInt { default: None, min: 1, max: 1000 },
        },
        ParamSpec {
            name: "nucmer_breaklen".into(),
            label: "Breakpoint distance".into(),
            help: "Advanced: how far the aligner searches for a continuation when an alignment breaks. Leave empty for the tool default. Changing this triggers a full re-alignment.".into(),
            layer: ParamLayer::Align,
            kind: ParamKind::OptionalInt { default: None, min: 1, max: 5000 },
        },
        ParamSpec {
            name: "dnadiff".into(),
            label: "Overall alignment report".into(),
            help: "Compute the summary statistics of the whole genome comparison.".into(),
            layer: ParamLayer::Postprocess,
            kind: ParamKind::Bool { default: true },
        },
    ]
}

#[derive(Debug, Clone, Serialize)]
pub struct ParamSpec {
    pub name: String,
    pub label: String,
    pub help: String,
    pub layer: ParamLayer,
    pub kind: ParamKind,
}

#[derive(Debug, Clone, Serialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum ParamKind {
    Int {
        default: i64,
        min: i64,
        max: i64,
    },
    Float {
        default: f64,
        min: f64,
        max: f64,
    },
    OptionalInt {
        default: Option<i64>,
        min: i64,
        max: i64,
    },
    Bool {
        default: bool,
    },
}

impl RunParams {
    /// Which layer a parameter belongs to (for cheap re-run detection).
    pub fn changed_params(&self, other: &RunParams) -> Vec<String> {
        let mut changed = Vec::new();
        if self.min_gap != other.min_gap {
            changed.push("min_gap".into());
        }
        if self.present_cov != other.present_cov {
            changed.push("present_cov".into());
        }
        if self.partial_cov != other.partial_cov {
            changed.push("partial_cov".into());
        }
        if self.blast_cov != other.blast_cov {
            changed.push("blast_cov".into());
        }
        if self.blast_pid != other.blast_pid {
            changed.push("blast_pid".into());
        }
        if self.blast_evalue != other.blast_evalue {
            changed.push("blast_evalue".into());
        }
        if self.dnadiff != other.dnadiff {
            changed.push("dnadiff".into());
        }
        changed
    }

    pub fn changed_align_params(&self, other: &RunParams) -> Vec<String> {
        let mut changed = Vec::new();
        if self.nucmer_minmatch != other.nucmer_minmatch {
            changed.push("nucmer_minmatch".into());
        }
        if self.nucmer_breaklen != other.nucmer_breaklen {
            changed.push("nucmer_breaklen".into());
        }
        changed
    }

    /// Stable string used in cache keys: only inputs and align-layer
    /// parameters invalidate a cached alignment.
    pub fn align_signature(&self) -> String {
        format!("l={:?};b={:?}", self.nucmer_minmatch, self.nucmer_breaklen)
    }
}

/// Validate a parameter map coming from JSON, returning field-level errors.
pub fn validate_params(params: &RunParams) -> Result<(), BTreeMap<String, String>> {
    use garde::Validate;
    let mut errors = BTreeMap::new();
    let ctx = ();
    match params.validate_with(&ctx) {
        Ok(()) => {}
        Err(e) => {
            for (path, err) in e.iter() {
                let field = path.to_string();
                if field.is_empty() {
                    continue;
                }
                errors.insert(field.clone(), friendly_range_error(&field, err.message()));
            }
        }
    }
    if params.partial_cov > params.present_cov {
        errors.insert(
            "partial_cov".into(),
            "The partial coverage threshold cannot be above the present threshold.".into(),
        );
    }
    if errors.is_empty() {
        Ok(())
    } else {
        Err(errors)
    }
}

fn friendly_range_error(field: &str, _raw: &str) -> String {
    let limit = match field {
        "min_gap" => "between 0 and 100000",
        "present_cov" | "partial_cov" | "blast_cov" | "blast_pid" => "between 0 and 100",
        "blast_evalue" => "between 1e-50 and 10",
        _ => "within the allowed range",
    };
    format!("Please enter a value {limit}.")
}
