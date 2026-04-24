//! Report generation for `MaxAlign` results.

use crate::alignment::AlignmentMetrics;
use crate::error::{Error, Result};
use crate::fasta::get_record_accession_string;
use crate::heuristic::HeuristicMethod;
use itertools::Itertools;
use serde::Serialize;
use std::collections::HashSet;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

#[derive(Serialize)]
struct JsonReport {
    run_options: RunOptions,
    statistics: Statistics,
    heuristic_iterations: Vec<IterationRecord>,
    refinement: Option<RefinementReport>,
    excluded_sequences: Vec<String>,
}

#[derive(Serialize)]
struct RunOptions {
    input_file: String,
    output_file: String,
    heuristic_method: u8,
    max_iterations: Option<u32>,
    improvement_threshold: f64,
    excluded_sequences_threshold: f64,
    keep_sequences: Vec<String>,
    retained_sequences_file: Option<String>,
    excluded_sequences_file: Option<String>,
}

#[derive(Serialize)]
struct Statistics {
    sequence_count: MetricChange,
    alignment_area: MetricChange,
    gap_free_columns: MetricChange,
    alignment_length: MetricChange,
}

#[derive(Serialize)]
struct MetricChange {
    before: usize,
    after: usize,
    change: i64,
}

#[derive(Serialize)]
struct IterationRecord {
    number: usize,
    excluded_this_round: usize,
    total_excluded: usize,
    gap_free_columns: usize,
    alignment_area: usize,
    excluded_sequences: Vec<String>,
}

#[derive(Serialize)]
struct RefinementReport {
    outcome: RefinementOutcome,
    bounded_by_excluded_sequences_threshold: bool,
    heuristic_metrics: AlignmentMetricsReport,
    final_metrics: AlignmentMetricsReport,
}

#[derive(Serialize)]
#[serde(rename_all = "snake_case")]
enum RefinementOutcome {
    HeuristicOptimal,
    HeuristicOptimalWithinExcludedSequencesThreshold,
    RetainedSequenceImprovement,
    AlignmentAreaImprovement,
}

#[derive(Serialize)]
struct AlignmentMetricsReport {
    sequence_count: usize,
    gap_free_columns: usize,
    alignment_area: usize,
    alignment_length: usize,
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
    pub retained_sequences: Option<PathBuf>,
    pub excluded_sequences: Option<PathBuf>,
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

/// Writes a detailed JSON report of `MaxAlign` results.
pub fn write_report(
    path: impl AsRef<Path>,
    config: &ReportConfig<'_>,
    data: &ReportData<'_>,
) -> Result<()> {
    let path = path.as_ref();
    let report = build_report(config, data);
    let serialized = serde_json::to_vec_pretty(&report).map_err(|e| Error::ReportSerialize {
        path: path.to_path_buf(),
        source: e,
    })?;

    let file = std::fs::File::create(path).map_err(|e| Error::ReportWrite {
        path: path.to_path_buf(),
        source: e,
    })?;
    let mut writer = BufWriter::new(file);
    writer
        .write_all(&serialized)
        .and_then(|()| writer.write_all(b"\n"))
        .and_then(|()| writer.flush())
        .map_err(|e| Error::ReportWrite {
            path: path.to_path_buf(),
            source: e,
        })
}

fn build_report(config: &ReportConfig<'_>, data: &ReportData<'_>) -> JsonReport {
    JsonReport {
        run_options: RunOptions {
            input_file: config.input_path.clone(),
            output_file: config.output_path.clone(),
            heuristic_method: config.heuristic_method as u8,
            max_iterations: (config.max_iterations != u32::MAX).then_some(config.max_iterations),
            improvement_threshold: config.improvement_threshold,
            excluded_sequences_threshold: config.excluded_seqs_threshold,
            keep_sequences: config.keep_sequence.to_vec(),
            retained_sequences_file: optional_path(config.retained_sequences.as_deref()),
            excluded_sequences_file: optional_path(config.excluded_sequences.as_deref()),
        },
        statistics: Statistics {
            sequence_count: metric_change(
                data.initial_metrics.sequence_count,
                data.final_metrics.sequence_count,
            ),
            alignment_area: metric_change(
                data.initial_metrics.alignment_area,
                data.final_metrics.alignment_area,
            ),
            gap_free_columns: metric_change(
                data.initial_metrics.gap_free_columns,
                data.final_metrics.gap_free_columns,
            ),
            alignment_length: metric_change(
                data.initial_metrics.alignment_length,
                data.final_metrics.alignment_length,
            ),
        },
        heuristic_iterations: iteration_records(
            data.iteration_data,
            data.initial_metrics,
            data.headers,
        ),
        refinement: refinement_report(config, data.heuristic_metrics, data.final_metrics),
        excluded_sequences: excluded_accessions(data.headers, data.excluded),
    }
}

fn optional_path(path: Option<&Path>) -> Option<String> {
    path.map(|path| path.display().to_string())
}

#[allow(clippy::cast_possible_wrap)]
fn metric_change(before: usize, after: usize) -> MetricChange {
    MetricChange {
        before,
        after,
        change: after as i64 - before as i64,
    }
}

fn iteration_records(
    iteration_data: &[(Vec<usize>, usize)],
    initial_metrics: &AlignmentMetrics,
    headers: &[Vec<u8>],
) -> Vec<IterationRecord> {
    let mut cumulative_excluded = 0;
    let mut records = Vec::with_capacity(iteration_data.len());

    for (i, (excluded_seqs, alignment_area)) in iteration_data.iter().enumerate() {
        cumulative_excluded += excluded_seqs.len();
        let remaining_seqs = initial_metrics.sequence_count - cumulative_excluded;
        let gap_free_columns = alignment_area.checked_div(remaining_seqs).unwrap_or(0);
        let excluded_sequences = excluded_seqs
            .iter()
            .map(|&idx| accession(headers, idx))
            .collect();

        records.push(IterationRecord {
            number: i + 1,
            excluded_this_round: excluded_seqs.len(),
            total_excluded: cumulative_excluded,
            gap_free_columns,
            alignment_area: *alignment_area,
            excluded_sequences,
        });
    }

    records
}

fn refinement_report(
    config: &ReportConfig<'_>,
    heuristic_metrics: &AlignmentMetrics,
    final_metrics: &AlignmentMetrics,
) -> Option<RefinementReport> {
    if !config.refinement {
        return None;
    }

    let bounded_by_excluded_sequences_threshold = config.excluded_seqs_threshold < 1.0;
    let outcome = if heuristic_metrics.alignment_area == final_metrics.alignment_area
        && final_metrics.sequence_count > heuristic_metrics.sequence_count
    {
        RefinementOutcome::RetainedSequenceImprovement
    } else if heuristic_metrics.alignment_area == final_metrics.alignment_area
        && bounded_by_excluded_sequences_threshold
    {
        RefinementOutcome::HeuristicOptimalWithinExcludedSequencesThreshold
    } else if heuristic_metrics.alignment_area == final_metrics.alignment_area {
        RefinementOutcome::HeuristicOptimal
    } else {
        RefinementOutcome::AlignmentAreaImprovement
    };

    Some(RefinementReport {
        outcome,
        bounded_by_excluded_sequences_threshold,
        heuristic_metrics: metrics_report(heuristic_metrics),
        final_metrics: metrics_report(final_metrics),
    })
}

fn metrics_report(metrics: &AlignmentMetrics) -> AlignmentMetricsReport {
    AlignmentMetricsReport {
        sequence_count: metrics.sequence_count,
        gap_free_columns: metrics.gap_free_columns,
        alignment_area: metrics.alignment_area,
        alignment_length: metrics.alignment_length,
    }
}

fn excluded_accessions(headers: &[Vec<u8>], excluded: &HashSet<usize>) -> Vec<String> {
    excluded
        .iter()
        .sorted_unstable()
        .map(|&idx| accession(headers, idx))
        .collect()
}

fn accession(headers: &[Vec<u8>], idx: usize) -> String {
    headers
        .get(idx)
        .and_then(|header| get_record_accession_string(header))
        .unwrap_or_default()
}
