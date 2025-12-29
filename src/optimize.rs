//! Branch-and-bound optimization algorithm.

use crate::alignment::AlignmentMetrics;
use crate::bitops::{
    bitwise_or, bitwise_or_assign, count_bits, count_bits_union, get_set_bit_indices,
};
use crate::heuristic::create_working_sets;
use log::debug;
use std::collections::HashSet;

const UNDECIDED: u8 = b'X';
const EXCLUDED: u8 = b'1';
const NOT_EXCLUDED: u8 = b'0';

/// Result of the branch-and-bound optimization.
pub struct BranchAndBoundResult {
    pub metrics: AlignmentMetrics,
    pub excluded: HashSet<usize>,
}

/// Runs the branch-and-bound algorithm to find the optimal solution.
#[must_use]
pub fn run_branch_and_bound(
    orig_sets: &[Vec<u8>],
    orig_gaps: &[Vec<u8>],
    metrics: &AlignmentMetrics,
    keep_pattern: &[bool],
    num_sequences: usize,
) -> BranchAndBoundResult {
    let kept_gaps = keep_pattern.iter().filter(|&&b| b).count();
    let gap_free_columns = metrics.alignment_length - orig_sets.len() - kept_gaps;

    let bb_state_excluded = HashSet::new();
    let bb_state_translation: Vec<usize> = (0..num_sequences).collect();

    let (working_sets, working_gaps) = create_working_sets(
        orig_sets,
        orig_gaps,
        &bb_state_excluded,
        &bb_state_translation,
        num_sequences,
    );

    let mut current_sets = working_sets;
    let mut current_gaps = working_gaps;

    crate::alignment::congruent_set_joining(
        &mut current_sets,
        &mut current_gaps,
        metrics.alignment_area,
        num_sequences,
        metrics.alignment_length,
    );

    crate::alignment::subset_joining(&current_sets, &mut current_gaps);

    let gap_columns = crate::alignment::set_elimination(
        &mut current_sets,
        &mut current_gaps,
        metrics.alignment_area,
        num_sequences,
        metrics.alignment_length,
        gap_free_columns,
    );

    let dislikes = find_dislikes(
        &current_sets,
        metrics.alignment_area,
        num_sequences,
        metrics.alignment_length,
        gap_columns,
    );

    let (ordered_sets, ordered_gaps, ordered_dislikes) =
        reorder_sets_for_search(&current_sets, &current_gaps, &dislikes);

    let (best_area, solutions) = branch_and_bound_search(
        &ordered_sets,
        &ordered_gaps,
        &ordered_dislikes,
        metrics.alignment_area,
        gap_free_columns,
        num_sequences,
    );

    extract_best_solution(solutions, best_area, num_sequences, metrics)
}

/// Performs the actual branch-and-bound search.
fn branch_and_bound_search(
    ordered_sets: &[Vec<u8>],
    ordered_gaps: &[Vec<u8>],
    ordered_dislikes: &[Vec<usize>],
    initial_best_area: usize,
    gap_free_columns: usize,
    num_sequences: usize,
) -> (usize, Vec<Vec<u8>>) {
    let sets_count = ordered_sets.len();
    let gap_vec_len = ordered_gaps.first().map_or(1, Vec::len);

    let mut suffix_unions = vec![vec![0u8; gap_vec_len]; sets_count + 1];
    for i in (0..sets_count).rev() {
        suffix_unions[i] = suffix_unions[i + 1].clone();
        bitwise_or_assign(&mut suffix_unions[i], &ordered_gaps[i]);
    }

    let mut stack = vec![(
        vec![UNDECIDED; sets_count],
        0usize,
        vec![0u8; num_sequences.div_ceil(8)],
        vec![0u8; gap_vec_len],
    )];

    let mut solutions = Vec::new();
    let mut best_area = initial_best_area;

    let union_sets_count_bits = |union_sets: &[u8]| count_bits(union_sets);

    while let Some((mut decisions, mut pointer, mut union_sets, mut union_gaps)) = stack.pop() {
        loop {
            let test_union_gaps_count = count_bits_union(&union_gaps, &suffix_unions[pointer]);
            let test_score = (gap_free_columns + test_union_gaps_count)
                * (num_sequences - union_sets_count_bits(&union_sets));

            if test_score < best_area {
                break;
            }

            while pointer < sets_count && decisions[pointer] != UNDECIDED {
                pointer += 1;
            }

            if pointer < sets_count {
                let set = &ordered_sets[pointer];
                let union_and_set = bitwise_or(&union_sets, set);

                if union_and_set == union_sets {
                    bitwise_or_assign(&mut union_gaps, &ordered_gaps[pointer]);
                    decisions[pointer] = EXCLUDED;
                    pointer += 1;
                    continue;
                }

                let mut decisions_not_excluded = decisions.clone();
                decisions_not_excluded[pointer] = NOT_EXCLUDED;
                stack.push((
                    decisions_not_excluded,
                    pointer + 1,
                    union_sets.clone(),
                    union_gaps.clone(),
                ));

                decisions[pointer] = EXCLUDED;
                union_sets = union_and_set;
                bitwise_or_assign(&mut union_gaps, &ordered_gaps[pointer]);

                for &bad in &ordered_dislikes[pointer] {
                    if bad > pointer {
                        decisions[bad] = NOT_EXCLUDED;
                    }
                }
                pointer += 1;
                continue;
            }

            let score = (gap_free_columns + count_bits(&union_gaps))
                * (num_sequences - count_bits(&union_sets));
            if score > best_area {
                best_area = score;
                solutions = vec![union_sets.clone()];
                debug!(
                    "Refinement algorithm improved the alignment: the area increased to {} with {} sequences",
                    best_area,
                    num_sequences - count_bits(&union_sets)
                );
            } else if score == best_area {
                solutions.push(union_sets.clone());
            }
            break;
        }
    }

    (best_area, solutions)
}

fn extract_best_solution(
    solutions: Vec<Vec<u8>>,
    best_area: usize,
    num_sequences: usize,
    metrics: &AlignmentMetrics,
) -> BranchAndBoundResult {
    if let Some(best_solution) = solutions.into_iter().min_by_key(|s| count_bits(s)) {
        let excluded_indices = get_set_bit_indices(&best_solution, num_sequences);
        let excluded: HashSet<usize> = excluded_indices.into_iter().collect();

        let remaining_seqs = num_sequences - excluded.len();
        if remaining_seqs > 0 {
            return BranchAndBoundResult {
                metrics: AlignmentMetrics::new(
                    remaining_seqs,
                    best_area / remaining_seqs,
                    best_area,
                    metrics.alignment_length,
                ),
                excluded,
            };
        }
    }

    BranchAndBoundResult {
        metrics: metrics.clone(),
        excluded: HashSet::new(),
    }
}

/// Finds pairs of sets that "dislike" each other (one is subset of other, or
/// their union would be too large to improve alignment area).
fn find_dislikes(
    sets: &[Vec<u8>],
    alignment_area: usize,
    sequence_count: usize,
    alignment_length: usize,
    gap_columns: usize,
) -> Vec<Vec<usize>> {
    let mut dislikes = vec![Vec::new(); sets.len()];
    let sets_count = sets.len();
    for i in 0..sets_count {
        let set_i = &sets[i];
        let set_i_bits = count_bits(set_i);
        for j in i + 1..sets_count {
            let set_j = &sets[j];
            let union_size = count_bits_union(set_i, set_j);
            if union_size == set_i_bits || union_size == count_bits(set_j) {
                dislikes[i].push(j);
                dislikes[j].push(i);
                continue;
            }

            if alignment_area > (alignment_length - gap_columns) * (sequence_count - union_size) {
                dislikes[i].push(j);
                dislikes[j].push(i);
            }
        }
    }
    dislikes
}

#[allow(clippy::type_complexity)]
fn reorder_sets_for_search(
    sets: &[Vec<u8>],
    gaps: &[Vec<u8>],
    dislikes: &[Vec<usize>],
) -> (Vec<Vec<u8>>, Vec<Vec<u8>>, Vec<Vec<usize>>) {
    let mut indices: Vec<usize> = (0..sets.len()).collect();

    indices.sort_by(|&a, &b| {
        dislikes[b]
            .len()
            .cmp(&dislikes[a].len())
            .then_with(|| count_bits(&sets[b]).cmp(&count_bits(&sets[a])))
            .then_with(|| count_bits(&gaps[b]).cmp(&count_bits(&gaps[a])))
    });

    let ordered_sets: Vec<Vec<u8>> = indices.iter().map(|&idx| sets[idx].clone()).collect();
    let ordered_gaps: Vec<Vec<u8>> = indices.iter().map(|&idx| gaps[idx].clone()).collect();

    let ordered_dislikes: Vec<Vec<usize>> = indices
        .iter()
        .map(|&idx| {
            dislikes[idx]
                .iter()
                .filter_map(|&d| indices.iter().position(|&x| x == d))
                .collect()
        })
        .collect();

    (ordered_sets, ordered_gaps, ordered_dislikes)
}
