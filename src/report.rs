//! Report generation for `MaxAlign` results.

use crate::alignment::AlignmentMetrics;
use crate::error::{Error, Result};
use crate::fasta::get_record_accession_string;
use crate::heuristic::HeuristicMethod;
use itertools::Itertools;
use markdown_tables::{MarkdownTableRow, as_table};
use std::collections::HashSet;
use std::io::{BufWriter, Write};
use std::path::Path;

struct RunOption {
    option: String,
    value: String,
}

impl MarkdownTableRow for RunOption {
    fn column_names() -> Vec<&'static str> {
        vec!["Option", "Value"]
    }

    fn column_values(&self) -> Vec<String> {
        vec![self.option.clone(), self.value.clone()]
    }
}

struct Statistic {
    metric: String,
    before: usize,
    after: usize,
    change: i64,
}

impl MarkdownTableRow for Statistic {
    fn column_names() -> Vec<&'static str> {
        vec!["Metric", "Before", "After", "Change"]
    }

    fn column_values(&self) -> Vec<String> {
        vec![
            self.metric.clone(),
            self.before.to_string(),
            self.after.to_string(),
            format!("{:+}", self.change),
        ]
    }
}

struct IterationRecord {
    number: usize,
    excluded_this_round: usize,
    total_excluded: usize,
    ungapped_columns: usize,
    alignment_area: usize,
}

impl MarkdownTableRow for IterationRecord {
    fn column_names() -> Vec<&'static str> {
        vec![
            "Iteration",
            "Excluded in this iteration",
            "Total excluded",
            "Ungapped columns",
            "Alignment area",
        ]
    }

    fn column_values(&self) -> Vec<String> {
        vec![
            self.number.to_string(),
            self.excluded_this_round.to_string(),
            self.total_excluded.to_string(),
            self.ungapped_columns.to_string(),
            self.alignment_area.to_string(),
        ]
    }
}

/// Configuration for generating a report.
#[derive(Debug)]
pub struct ReportConfig<'a> {
    pub input_path: String,
    pub output_path: String,
    pub heuristic_method: HeuristicMethod,
    pub max_iterations: u32,
    pub improvement_threshold: f64,
    pub excluded_seqs_threshold: f64,
    pub refinement: bool,
    pub keep_sequence: &'a [String],
    pub retained_sequences: Option<String>,
    pub excluded_sequences: Option<String>,
}

/// Data for generating a report.
pub struct ReportData<'a> {
    pub initial_metrics: &'a AlignmentMetrics,
    pub heuristic_metrics: &'a AlignmentMetrics,
    pub final_metrics: &'a AlignmentMetrics,
    pub iteration_data: &'a [(Vec<usize>, usize)],
    pub headers: &'a [Vec<u8>],
    pub excluded: &'a HashSet<usize>,
}

/// Writes a detailed report of `MaxAlign` results.
#[allow(clippy::cast_possible_wrap)]
pub fn write_report(
    path: impl AsRef<Path>,
    config: &ReportConfig<'_>,
    data: &ReportData<'_>,
) -> Result<()> {
    let path = path.as_ref();
    let file = std::fs::File::create(path).map_err(|e| Error::ReportWrite {
        path: path.to_path_buf(),
        source: e,
    })?;
    let mut writer = BufWriter::new(file);

    write_header(&mut writer, path)?;
    write_options_section(&mut writer, config, path)?;
    write_statistics_section(&mut writer, data.initial_metrics, data.final_metrics, path)?;
    write_iterations_section(&mut writer, data.iteration_data, data.initial_metrics, path)?;
    write_refinement_section(
        &mut writer,
        config,
        data.heuristic_metrics,
        data.final_metrics,
        path,
    )?;
    write_excluded_section(&mut writer, data.headers, data.excluded, path)?;

    writer.flush().map_err(|e| Error::ReportWrite {
        path: path.to_path_buf(),
        source: e,
    })?;

    Ok(())
}

macro_rules! write_err {
    ($path:expr) => {
        |e| Error::ReportWrite {
            path: $path.to_path_buf(),
            source: e,
        }
    };
}

fn write_header(writer: &mut impl Write, path: &Path) -> Result<()> {
    writeln!(writer, "# MaxAlign Results\n").map_err(write_err!(path))
}

fn write_options_section(
    writer: &mut impl Write,
    config: &ReportConfig<'_>,
    report_path: &Path,
) -> Result<()> {
    writeln!(writer, "## Run options\n").map_err(write_err!(report_path))?;

    let max_iter_str = if config.max_iterations == u32::MAX {
        "unlimited".to_string()
    } else {
        config.max_iterations.to_string()
    };

    let mut options = vec![
        RunOption {
            option: "Input file".to_string(),
            value: config.input_path.clone(),
        },
        RunOption {
            option: "Output file".to_string(),
            value: config.output_path.clone(),
        },
        RunOption {
            option: "Heuristic method".to_string(),
            value: config.heuristic_method.to_string(),
        },
        RunOption {
            option: "Max iterations".to_string(),
            value: max_iter_str,
        },
        RunOption {
            option: "Improvement threshold".to_string(),
            value: config.improvement_threshold.to_string(),
        },
        RunOption {
            option: "Excluded sequences threshold".to_string(),
            value: config.excluded_seqs_threshold.to_string(),
        },
        RunOption {
            option: "Refinement".to_string(),
            value: config.refinement.to_string(),
        },
        RunOption {
            option: "Keep sequences".to_string(),
            value: if config.keep_sequence.is_empty() {
                String::new()
            } else {
                config.keep_sequence.join(", ")
            },
        },
    ];

    if let Some(ref retained) = config.retained_sequences {
        options.push(RunOption {
            option: "Retained sequences file".to_string(),
            value: retained.clone(),
        });
    }

    if let Some(ref excluded) = config.excluded_sequences {
        options.push(RunOption {
            option: "Excluded sequences file".to_string(),
            value: excluded.clone(),
        });
    }

    options.push(RunOption {
        option: "Report file".to_string(),
        value: report_path.display().to_string(),
    });

    writeln!(writer, "{}", as_table(&options)).map_err(write_err!(report_path))
}

#[allow(clippy::cast_possible_wrap)]
fn write_statistics_section(
    writer: &mut impl Write,
    initial_metrics: &AlignmentMetrics,
    final_metrics: &AlignmentMetrics,
    path: &Path,
) -> Result<()> {
    writeln!(writer, "## Statistics\n").map_err(write_err!(path))?;

    let sequences_change =
        final_metrics.sequence_count as i64 - initial_metrics.sequence_count as i64;
    let area_change = final_metrics.alignment_area as i64 - initial_metrics.alignment_area as i64;
    let freecols_change =
        final_metrics.gap_free_columns as i64 - initial_metrics.gap_free_columns as i64;
    let totalcols_change =
        final_metrics.alignment_length as i64 - initial_metrics.alignment_length as i64;

    let statistics = vec![
        Statistic {
            metric: "Number of sequences".to_string(),
            before: initial_metrics.sequence_count,
            after: final_metrics.sequence_count,
            change: sequences_change,
        },
        Statistic {
            metric: "Alignment area".to_string(),
            before: initial_metrics.alignment_area,
            after: final_metrics.alignment_area,
            change: area_change,
        },
        Statistic {
            metric: "Ungapped columns".to_string(),
            before: initial_metrics.gap_free_columns,
            after: final_metrics.gap_free_columns,
            change: freecols_change,
        },
        Statistic {
            metric: "Total columns".to_string(),
            before: initial_metrics.alignment_length,
            after: final_metrics.alignment_length,
            change: totalcols_change,
        },
    ];

    writeln!(writer, "{}", as_table(&statistics)).map_err(write_err!(path))
}

fn write_iterations_section(
    writer: &mut impl Write,
    iteration_data: &[(Vec<usize>, usize)],
    initial_metrics: &AlignmentMetrics,
    path: &Path,
) -> Result<()> {
    writeln!(writer, "## Heuristic iterations\n").map_err(write_err!(path))?;

    if iteration_data.is_empty() {
        writeln!(
            writer,
            "No iterations performed. Alignment could not be improved.\n"
        )
        .map_err(write_err!(path))
    } else {
        let mut cumulative_excluded = 0;
        let mut iterations = Vec::new();
        for (i, (excluded_seqs, align_area)) in iteration_data.iter().enumerate() {
            cumulative_excluded += excluded_seqs.len();
            let remaining_seqs = initial_metrics.sequence_count - cumulative_excluded;
            let freecols = if remaining_seqs > 0 {
                align_area / remaining_seqs
            } else {
                0
            };
            iterations.push(IterationRecord {
                number: i + 1,
                excluded_this_round: excluded_seqs.len(),
                total_excluded: cumulative_excluded,
                ungapped_columns: freecols,
                alignment_area: *align_area,
            });
        }
        writeln!(writer, "{}", as_table(&iterations)).map_err(write_err!(path))
    }
}

fn write_refinement_section(
    writer: &mut impl Write,
    config: &ReportConfig<'_>,
    heuristic_metrics: &AlignmentMetrics,
    final_metrics: &AlignmentMetrics,
    path: &Path,
) -> Result<()> {
    if !config.refinement {
        return Ok(());
    }

    writeln!(writer, "## Refinement\n").map_err(write_err!(path))?;

    if heuristic_metrics.alignment_area == final_metrics.alignment_area {
        writeln!(
            writer,
            "The solution found with the heuristic method is optimal, as \
             one determined by the branch-and-bound algorithm. The alignment \
             area remains {}.\n",
            heuristic_metrics.alignment_area
        )
        .map_err(write_err!(path))
    } else {
        writeln!(
            writer,
            "The heuristic solution was improved by the branch-and-bound algorithm. \
             The alignment area increased from {} to {}.\n",
            heuristic_metrics.alignment_area, final_metrics.alignment_area
        )
        .map_err(write_err!(path))
    }
}

fn write_excluded_section(
    writer: &mut impl Write,
    headers: &[Vec<u8>],
    excluded: &HashSet<usize>,
    path: &Path,
) -> Result<()> {
    writeln!(writer, "## Excluded sequences\n").map_err(write_err!(path))?;

    if excluded.is_empty() {
        writeln!(writer, "No sequences were excluded.").map_err(write_err!(path))
    } else {
        // Write excluded sequences as a simple bullet list (no indices).
        for name in excluded
            .iter()
            .sorted_unstable()
            .map(|&idx| get_record_accession_string(&headers[idx]).unwrap_or_default())
        {
            writeln!(writer, "- {}", name).map_err(write_err!(path))?;
        }
        Ok(())
    }
}
