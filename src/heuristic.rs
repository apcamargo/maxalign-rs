//! Heuristic algorithm for sequence exclusion.

use crate::alignment::{AlignmentMetrics, SetData, congruent_set_joining, subset_joining};
use crate::bitops::{
    bitwise_or, count_bits, count_bits_union, count_bits_union_triple, get_set_bit_indices,
    pack_bools_to_bits,
};
use log::info;
use std::collections::HashSet;

/// The heuristic method to use for finding sequences to exclude.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[allow(clippy::enum_variant_names)]
pub enum HeuristicMethod {
    NoSynergy = 1,
    #[default]
    PairwiseSynergy = 2,
    TripleSynergy = 3,
}

impl std::fmt::Display for HeuristicMethod {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", *self as u8)
    }
}

impl std::str::FromStr for HeuristicMethod {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "1" => Ok(Self::NoSynergy),
            "2" => Ok(Self::PairwiseSynergy),
            "3" => Ok(Self::TripleSynergy),
            _ => Err(format!(
                "invalid heuristic method '{s}': must be 1, 2, or 3"
            )),
        }
    }
}

/// Configuration for the heuristic algorithm.
#[derive(Debug, Clone)]
pub struct HeuristicConfig {
    pub method: HeuristicMethod,
    pub max_iterations: u32,
    pub improvement_threshold: f64,
    pub excluded_seqs_threshold: f64,
}

/// Runs the heuristic algorithm to find sequences to exclude.
#[allow(clippy::cast_precision_loss)]
pub fn run_heuristic(
    state: &mut SetData,
    metrics: &mut AlignmentMetrics,
    config: &HeuristicConfig,
    keep_pattern: &[bool],
    num_orig_seqs: usize,
) -> Vec<(Vec<usize>, usize)> {
    let kept_gaps_count = keep_pattern.iter().filter(|&&b| b).count();
    let mut iteration_data = Vec::new();
    let mut iterations_count: u32 = 0;

    loop {
        if iterations_count >= config.max_iterations {
            break;
        }

        let (working_sets, working_gaps) = create_working_sets(
            &state.sets,
            &state.gaps,
            &state.excluded,
            &state.translation,
            num_orig_seqs,
        );

        let sequence_count = state.translation.len();
        let gap_free_columns = metrics.alignment_length - working_sets.len() - kept_gaps_count;

        metrics.sequence_count = sequence_count;
        metrics.gap_free_columns = gap_free_columns;

        let mut current_sets = working_sets;
        let mut current_gaps = working_gaps;

        congruent_set_joining(
            &mut current_sets,
            &mut current_gaps,
            metrics.alignment_area,
            metrics.sequence_count,
            metrics.alignment_length,
        );

        subset_joining(&current_sets, &mut current_gaps);

        if current_sets.is_empty() {
            break;
        }

        let (best_set, new_alignment_area) = find_greatest_impact_set(
            &current_sets,
            &current_gaps,
            metrics.alignment_area,
            sequence_count,
            gap_free_columns,
            config.method,
        );

        if config.improvement_threshold != 0.0 && metrics.alignment_area != 0 {
            let improvement = (new_alignment_area as f64 - metrics.alignment_area as f64)
                / metrics.alignment_area as f64;
            if improvement < config.improvement_threshold {
                info!(
                    "Early stopping: relative improvement ({:.4}) is below threshold ({:.4})",
                    improvement, config.improvement_threshold
                );
                break;
            }
        }

        let excluded_fraction = state.excluded.len() as f64 / num_orig_seqs as f64;
        if excluded_fraction >= config.excluded_seqs_threshold {
            info!(
                "Early stopping: excluded sequence fraction ({:.4}) reached threshold ({:.4})",
                excluded_fraction, config.excluded_seqs_threshold
            );
            break;
        }

        if metrics.alignment_area >= new_alignment_area {
            break;
        }

        let excluded_indices = get_set_bit_indices(&best_set, sequence_count);
        let mut exseq = Vec::with_capacity(excluded_indices.len());
        for pointer in excluded_indices {
            let orig_idx = state.translation[pointer];
            state.excluded.insert(orig_idx);
            exseq.push(orig_idx);
        }

        state.translation = (0..num_orig_seqs)
            .filter(|&i| !state.excluded.contains(&i))
            .collect();

        metrics.sequence_count = state.translation.len();
        metrics.alignment_area = new_alignment_area;
        iteration_data.push((exseq, new_alignment_area));
        iterations_count += 1;
    }

    iteration_data
}

/// Finds the set that, when excluded, provides the greatest improvement.
#[allow(clippy::cast_precision_loss, clippy::float_cmp)]
fn find_greatest_impact_set(
    sets: &[Vec<u8>],
    gaps: &[Vec<u8>],
    current_area: usize,
    sequence_count: usize,
    gap_free_columns: usize,
    method: HeuristicMethod,
) -> (Vec<u8>, usize) {
    let mut best_set = Vec::new();
    let mut best_impact = 0;
    let mut best_efficiency = -1.0;
    let mut best_gap_count = 0;

    let mut evaluate_candidate =
        |set_size: usize, gap_count: usize, candidate_fn: &dyn Fn() -> Vec<u8>| {
            let this_impact = (sequence_count - set_size) * (gap_free_columns + gap_count);
            let this_efficiency = (this_impact as f64 - current_area as f64) / set_size as f64;

            if this_efficiency > best_efficiency
                || (this_efficiency == best_efficiency && gap_count >= best_gap_count)
            {
                best_efficiency = this_efficiency;
                best_impact = this_impact;
                best_set = candidate_fn();
                best_gap_count = gap_count;
            }
        };

    for (i, set_i) in sets.iter().enumerate() {
        evaluate_candidate(count_bits(set_i), count_bits(&gaps[i]), &|| set_i.to_vec());

        if method as u8 >= 2 {
            for j in 0..i {
                evaluate_candidate(
                    count_bits_union(set_i, &sets[j]),
                    count_bits_union(&gaps[i], &gaps[j]),
                    &|| bitwise_or(set_i, &sets[j]),
                );

                if method as u8 >= 3 {
                    for k in 0..j {
                        evaluate_candidate(
                            count_bits_union_triple(set_i, &sets[j], &sets[k]),
                            count_bits_union_triple(&gaps[i], &gaps[j], &gaps[k]),
                            &|| bitwise_or(&bitwise_or(set_i, &sets[j]), &sets[k]),
                        );
                    }
                }
            }
        }
    }

    (best_set, best_impact)
}

/// Creates working sets by filtering out excluded sequences.
pub fn create_working_sets(
    orig_sets: &[Vec<u8>],
    orig_gaps: &[Vec<u8>],
    excluded: &HashSet<usize>,
    translation: &[usize],
    num_orig_seqs: usize,
) -> (Vec<Vec<u8>>, Vec<Vec<u8>>) {
    let mut working_sets = Vec::new();
    let mut working_gaps = Vec::new();

    let included: Vec<bool> = (0..num_orig_seqs)
        .map(|idx| !excluded.contains(&idx))
        .collect();

    for (i, orig_set) in orig_sets.iter().enumerate() {
        let mut bools = Vec::with_capacity(translation.len());
        let mut has_any_gap = false;

        for (idx, &is_included) in included.iter().enumerate() {
            if is_included {
                let byte_idx = idx / 8;
                let bit_idx = idx % 8;
                let is_gap = byte_idx < orig_set.len() && (orig_set[byte_idx] >> bit_idx) & 1 == 1;
                bools.push(is_gap);
                has_any_gap |= is_gap;
            }
        }

        if has_any_gap {
            working_sets.push(pack_bools_to_bits(&bools));
            working_gaps.push(orig_gaps[i].clone());
        }
    }

    (working_sets, working_gaps)
}
