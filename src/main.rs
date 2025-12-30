mod alignment;
mod bitops;
mod error;
mod fasta;
mod heuristic;
mod optimize;
mod output;
mod report;

use crate::alignment::{
    AlignmentMetrics, SetData, create_gap_matrix, create_sets, remove_all_gap_columns,
};
use crate::error::{Error, Result};
use crate::fasta::parse_fasta;
use crate::heuristic::{HeuristicConfig, HeuristicMethod, run_heuristic};
use crate::optimize::run_branch_and_bound;
use crate::output::{write_fasta, write_headers_list};
use crate::report::{ReportConfig, ReportData, write_report};
use clap::{
    CommandFactory, Parser,
    builder::styling::{AnsiColor, Style, Styles},
};
use clio::{Input, Output};
use env_logger::Builder;
use itertools::Itertools;
use log::{LevelFilter, debug, info};
use std::io::{IsTerminal, Write};
use std::process::ExitCode;

const STYLES: Styles = Styles::styled()
    .header(AnsiColor::Cyan.on_default().bold())
    .usage(AnsiColor::Yellow.on_default().bold())
    .literal(AnsiColor::Yellow.on_default().bold())
    .placeholder(Style::new().dimmed());

fn parse_max_iterations(s: &str) -> std::result::Result<u32, String> {
    if s == "-1" {
        return Ok(u32::MAX);
    }
    s.parse()
        .map_err(|_| format!("`{s}` is not a valid iteration count"))
}

fn parse_threshold(s: &str) -> std::result::Result<f64, String> {
    let v: f64 = s.parse().map_err(|_| format!("`{s}` isn't a number"))?;
    if v < 0.0 {
        Err("value must be non-negative".to_string())
    } else {
        Ok(v)
    }
}

#[derive(Parser)]
#[command(version, about, styles = STYLES, max_term_width = 88)]
struct Cli {
    /// Input FASTA file
    #[arg(default_value = "-")]
    input: Input,

    /// Output FASTA file
    #[arg(default_value = "-")]
    output: Output,

    /// Heuristic method: 1 (no synergy), 2 (pairwise synergy), 3 (three-way synergy)
    #[arg(short = 'm', long, default_value = "2", value_parser = clap::value_parser!(HeuristicMethod))]
    heuristic_method: HeuristicMethod,

    /// Maximum number of iterations (-1 for unlimited iterations)
    #[arg(short = 'i', long, default_value = "-1", value_parser = parse_max_iterations)]
    max_iterations: u32,

    /// Perform refinement using the branch-and-bound algorithm to find the optimal solution
    #[arg(short = 'o', long, default_value = "false")]
    refinement: bool,

    /// Stop iterating if the relative improvement is below this threshold
    #[arg(short = 't', long, default_value = "0.0", value_parser = parse_threshold)]
    improvement_threshold: f64,

    /// Stop iterating if the fraction of excluded sequences is above this threshold
    #[arg(short = 's', long, default_value = "1.0", value_parser = parse_threshold)]
    excluded_seqs_threshold: f64,

    /// Sequence to always retain (can be specified multiple times)
    #[arg(short = 'k', long)]
    keep_sequence: Vec<String>,

    /// Report file path
    #[arg(short = 'r', long)]
    report: Option<String>,

    /// Write a list of retained sequences to file
    #[arg(long)]
    retained_sequences: Option<String>,

    /// Write a list of excluded sequences to file
    #[arg(long)]
    excluded_sequences: Option<String>,

    /// Verbosity level (-v for normal logging, -vv for detailed logging)
    #[arg(short = 'v', long, action = clap::ArgAction::Count)]
    verbosity: u8,
}

fn setup_logging(verbosity: u8) {
    let level = match verbosity {
        0 => LevelFilter::Off,
        1 => LevelFilter::Info,
        2 => LevelFilter::Debug,
        _ => LevelFilter::Trace,
    };

    Builder::new()
        .filter_level(level)
        .format(|buf, record| writeln!(buf, "[{}] {}", buf.timestamp(), record.args()))
        .init();
}

#[allow(clippy::too_many_lines)]
fn run(cli: &Cli) -> Result<()> {
    if cli.input.is_std() && std::io::stdin().is_terminal() {
        #[allow(clippy::unwrap_used)]
        Cli::command().print_help().unwrap();
        return Ok(());
    }

    let sequence_data = match parse_fasta(&cli.input, &cli.keep_sequence) {
        Ok(data) => data,
        Err(Error::EmptyInput) if cli.input.is_std() => {
            #[allow(clippy::unwrap_used)]
            Cli::command().print_help().unwrap();
            return Ok(());
        }
        Err(e) => return Err(e),
    };

    let num_sequences = sequence_data.headers.len();

    let mut sequences = sequence_data.sequences;
    for seq in &mut sequences {
        seq.resize(sequence_data.longest_length, b'-');
    }

    let gap_matrix = create_gap_matrix(&sequences, sequence_data.longest_length);
    let (orig_sets, orig_gaps, keep_pattern) = create_sets(
        &gap_matrix,
        &sequence_data.keep_indices,
        sequence_data.longest_length,
    );

    let kept_gaps_count = keep_pattern.iter().filter(|&&b| b).count();
    let initial_gap_free_columns = sequence_data.longest_length - orig_sets.len() - kept_gaps_count;

    let initial_metrics = AlignmentMetrics::new(
        num_sequences,
        initial_gap_free_columns,
        initial_gap_free_columns * num_sequences,
        sequence_data.longest_length,
    );

    info!(
        "Loaded input alignment (sequences: {}, length: {}, initial area: {})",
        initial_metrics.sequence_count,
        initial_metrics.alignment_length,
        initial_metrics.alignment_area
    );

    let mut metrics = initial_metrics.clone();
    let mut state = SetData::new(orig_sets.clone(), orig_gaps.clone(), num_sequences);

    let heuristic_config = HeuristicConfig {
        method: cli.heuristic_method,
        max_iterations: cli.max_iterations,
        improvement_threshold: cli.improvement_threshold,
        excluded_seqs_threshold: cli.excluded_seqs_threshold,
    };

    info!(
        "Processing alignment (heuristic method: {})",
        cli.heuristic_method
    );
    let iteration_data = run_heuristic(
        &mut state,
        &mut metrics,
        &heuristic_config,
        &keep_pattern,
        num_sequences,
    );

    for (iter, (exseq, area)) in iteration_data.iter().enumerate() {
        let names = exseq
            .iter()
            .map(|&idx| {
                crate::fasta::get_record_accession_string(&sequence_data.headers[idx])
                    .unwrap_or_default()
            })
            .format(", ")
            .to_string();
        debug!(
            "Iteration {}: alignment area is {}, {} sequence(s) excluded ({})",
            iter + 1,
            area,
            exseq.len(),
            names
        );
    }

    let heuristic_metrics = metrics.clone();
    let mut final_excluded = state.excluded.clone();
    let mut final_metrics = metrics;

    if cli.refinement {
        info!(
            "Starting refinement using the branch-and-bound algorithm to find the optimal solution"
        );
        let bb_result = run_branch_and_bound(
            &orig_sets,
            &orig_gaps,
            &heuristic_metrics,
            &keep_pattern,
            num_sequences,
        );
        if bb_result.metrics.alignment_area > final_metrics.alignment_area {
            final_metrics = bb_result.metrics;
            final_excluded = bb_result.excluded;
        }
    }

    let excluded_count = initial_metrics.sequence_count - final_metrics.sequence_count;
    if excluded_count == 0 {
        info!(
            "No sequences were excluded. Alignment area remained {} ({} sequences)",
            initial_metrics.alignment_area, initial_metrics.sequence_count
        );
    } else {
        info!(
            "A total of {} sequences were excluded. Alignment area improved by {} (from {} to {})",
            excluded_count,
            final_metrics.alignment_area - initial_metrics.alignment_area,
            initial_metrics.alignment_area,
            final_metrics.alignment_area
        );
    }

    let (final_sequences, final_headers) =
        remove_all_gap_columns(&sequences, &sequence_data.headers, &final_excluded);

    if !final_sequences.is_empty()
        && let Some(final_alignment_length) = final_sequences.iter().map(Vec::len).max()
    {
        final_metrics.alignment_length = final_alignment_length;
        final_metrics.gap_free_columns =
            final_metrics.alignment_area / final_metrics.sequence_count;
    }

    let mut output = cli.output.clone();
    let output_name = if output.is_std() {
        "stdout".to_string()
    } else {
        output.path().to_string_lossy().into_owned()
    };

    write_fasta(&final_sequences, &final_headers, &mut output)?;
    info!("Output written to {}", output_name);

    if let Some(ref report_path) = cli.report {
        let input_path = if cli.input.is_std() {
            "<stdin>".to_string()
        } else {
            cli.input.path().to_string_lossy().to_string()
        };

        let output_path = if cli.output.is_std() {
            "<stdout>".to_string()
        } else {
            cli.output.path().to_string_lossy().to_string()
        };

        let config = ReportConfig {
            input_path,
            output_path,
            heuristic_method: cli.heuristic_method,
            max_iterations: cli.max_iterations,
            improvement_threshold: cli.improvement_threshold,
            excluded_seqs_threshold: cli.excluded_seqs_threshold,
            refinement: cli.refinement,
            keep_sequence: &cli.keep_sequence,
            retained_sequences: cli.retained_sequences.clone(),
            excluded_sequences: cli.excluded_sequences.clone(),
        };

        let data = ReportData {
            initial_metrics: &initial_metrics,
            heuristic_metrics: &heuristic_metrics,
            final_metrics: &final_metrics,
            iteration_data: &iteration_data,
            headers: &sequence_data.headers,
            excluded: &final_excluded,
        };

        write_report(report_path, &config, &data)?;
        info!("Report written to {}", report_path);
    }

    if let Some(ref path) = cli.retained_sequences {
        write_headers_list(path, &sequence_data.headers, &final_excluded, true)?;
        info!("List of retained sequences written to {}", path);
    }
    if let Some(ref path) = cli.excluded_sequences {
        write_headers_list(path, &sequence_data.headers, &final_excluded, false)?;
        info!("List of excluded sequences written to {}", path);
    }

    Ok(())
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    setup_logging(cli.verbosity);

    if let Err(e) = run(&cli) {
        eprintln!("Error: {e}");
        ExitCode::FAILURE
    } else {
        ExitCode::SUCCESS
    }
}
