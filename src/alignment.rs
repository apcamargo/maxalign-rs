//! Alignment data structures and operations.
//!
//! This module provides the core data structures for representing sequence
//! alignments and the operations needed to analyze gap patterns and compute
//! alignment metrics.

use crate::bitops::{bitwise_or_assign, count_bits, set_bit};
use std::collections::{HashMap, HashSet};

#[inline]
pub const fn is_gap_char(byte: u8) -> bool {
    byte == b'-' || byte == b'.'
}

/// Holds the current state of set data during optimization.
#[derive(Clone)]
pub struct SetData {
    pub sets: Vec<Vec<u8>>,
    pub gaps: Vec<Vec<u8>>,
    pub translation: Vec<usize>,
    pub excluded: HashSet<usize>,
}

impl SetData {
    #[must_use]
    pub fn new(sets: Vec<Vec<u8>>, gaps: Vec<Vec<u8>>, num_sequences: usize) -> Self {
        Self {
            sets,
            gaps,
            translation: (0..num_sequences).collect(),
            excluded: HashSet::new(),
        }
    }
}

/// Metrics describing the current state of an alignment.
#[derive(Clone, Debug, Default)]
pub struct AlignmentMetrics {
    pub sequence_count: usize,
    pub gap_free_columns: usize,
    pub alignment_area: usize,
    pub alignment_length: usize,
}

impl AlignmentMetrics {
    #[must_use]
    pub const fn new(
        sequence_count: usize,
        gap_free_columns: usize,
        alignment_area: usize,
        alignment_length: usize,
    ) -> Self {
        Self {
            sequence_count,
            gap_free_columns,
            alignment_area,
            alignment_length,
        }
    }
}

/// Creates a gap matrix from sequences where `gap_matrix[seq][col]` is `true`
/// if sequence `seq` has a gap at column `col`.
#[must_use]
pub fn create_gap_matrix(sequences: &[Vec<u8>], alignment_length: usize) -> Vec<Vec<bool>> {
    sequences
        .iter()
        .map(|seq| {
            let mut row = vec![true; alignment_length];
            for (col_idx, &byte) in seq.iter().enumerate() {
                if !is_gap_char(byte) {
                    row[col_idx] = false;
                }
            }
            row
        })
        .collect()
}

/// Creates gap pattern sets from a gap matrix, grouping columns by their gap pattern
/// and creating bit-packed representations for efficient manipulation.
#[must_use]
pub fn create_sets(
    gap_matrix: &[Vec<bool>],
    keep_indices: &HashSet<usize>,
    alignment_length: usize,
) -> (Vec<Vec<u8>>, Vec<Vec<u8>>, Vec<bool>) {
    let num_seqs = gap_matrix.len();
    let mut keep_pattern = vec![false; alignment_length];
    for &keep_seq_idx in keep_indices {
        for (col_idx, &is_gap) in gap_matrix[keep_seq_idx].iter().enumerate() {
            if is_gap {
                keep_pattern[col_idx] = true;
            }
        }
    }

    let bytes_per_col = num_seqs.div_ceil(8);
    let mut flat_sets = vec![0u8; alignment_length * bytes_per_col];
    let mut has_gap_col = vec![false; alignment_length];

    for (seq_idx, row) in gap_matrix.iter().enumerate() {
        let byte_offset = seq_idx / 8;
        let bit_mask = 1u8 << (seq_idx % 8);
        for (col_idx, &is_gap) in row.iter().enumerate() {
            if is_gap && !keep_pattern[col_idx] {
                flat_sets[col_idx * bytes_per_col + byte_offset] |= bit_mask;
                has_gap_col[col_idx] = true;
            }
        }
    }

    let mut sets = Vec::new();
    let mut gaps = Vec::new();

    for (col_idx, &has_gap) in has_gap_col.iter().enumerate() {
        if has_gap {
            let start = col_idx * bytes_per_col;
            let end = start + bytes_per_col;
            sets.push(flat_sets[start..end].to_vec());

            let mut gap_vec = vec![0u8; alignment_length.div_ceil(8)];
            set_bit(&mut gap_vec, col_idx);
            gaps.push(gap_vec);
        }
    }

    (sets, gaps, keep_pattern)
}

/// Joins congruent (identical) sets and removes sets that cannot improve the alignment.
/// Returns the number of gap columns that were removed.
pub fn congruent_set_joining(
    sets: &mut Vec<Vec<u8>>,
    gaps: &mut Vec<Vec<u8>>,
    alignment_area: usize,
    sequence_count: usize,
    alignment_length: usize,
) -> usize {
    let mut gap_columns = 0;
    let mut to_remove = HashSet::new();

    for (i, set) in sets.iter().enumerate() {
        let size_i = count_bits(set);
        if alignment_area > alignment_length * (sequence_count - size_i) {
            to_remove.insert(i);
            gap_columns += 1;
        }
    }

    let mut pattern_to_idx: HashMap<&[u8], usize> = HashMap::new();
    for i in (0..sets.len()).rev() {
        if to_remove.contains(&i) {
            continue;
        }

        if let Some(&last_idx) = pattern_to_idx.get(&sets[i] as &[u8]) {
            let gap_i = gaps[i].clone();
            bitwise_or_assign(&mut gaps[last_idx], &gap_i);
            to_remove.insert(i);
        } else {
            pattern_to_idx.insert(&sets[i], i);
        }
    }

    remove_indices_from_parallel_vecs(sets, gaps, to_remove);

    gap_columns
}

/// Propagates gap column benefits from subsets to their supersets.
pub fn subset_joining(sets: &[Vec<u8>], gaps: &mut [Vec<u8>]) {
    let mut merges: Vec<(usize, Vec<u8>)> = Vec::new();

    for i in 0..sets.len() {
        for j in 0..sets.len() {
            if i == j {
                continue;
            }
            if is_subset_of(&sets[j], &sets[i]) {
                merges.push((i, gaps[j].clone()));
            }
        }
    }

    for (target_idx, source_gap) in merges {
        bitwise_or_assign(&mut gaps[target_idx], &source_gap);
    }
}

#[inline]
fn is_subset_of(a: &[u8], b: &[u8]) -> bool {
    a.iter().zip(b).all(|(&x, &y)| (x & y) == x)
}

fn remove_indices_from_parallel_vecs(
    sets: &mut Vec<Vec<u8>>,
    gaps: &mut Vec<Vec<u8>>,
    to_remove: HashSet<usize>,
) {
    let mut indices: Vec<_> = to_remove.into_iter().collect();
    indices.sort_unstable_by(|a, b| b.cmp(a));
    for idx in indices {
        sets.swap_remove(idx);
        gaps.swap_remove(idx);
    }
}

/// Eliminates sets that cannot lead to an improvement in alignment area.
/// Returns the final number of gap columns.
pub fn set_elimination(
    sets: &mut Vec<Vec<u8>>,
    gaps: &mut Vec<Vec<u8>>,
    alignment_area: usize,
    sequence_count: usize,
    alignment_length: usize,
    gap_free_columns: usize,
) -> usize {
    let mut current_gap_columns = get_gap_columns(gaps, alignment_length, gap_free_columns);
    loop {
        let mut to_remove = HashSet::new();
        for (i, set) in sets.iter().enumerate() {
            let set_size = count_bits(set);
            if alignment_area
                > (alignment_length - current_gap_columns) * (sequence_count - set_size)
            {
                to_remove.insert(i);
            }
        }

        if to_remove.is_empty() {
            break;
        }

        remove_indices_from_parallel_vecs(sets, gaps, to_remove);

        let next_gap_columns = get_gap_columns(gaps, alignment_length, gap_free_columns);
        if next_gap_columns == current_gap_columns {
            break;
        }
        current_gap_columns = next_gap_columns;
    }
    current_gap_columns
}

/// Calculates the number of gap columns from gap indicators.
#[must_use]
pub fn get_gap_columns(
    gaps: &[Vec<u8>],
    alignment_length: usize,
    gap_free_columns: usize,
) -> usize {
    if gaps.is_empty() {
        return 0;
    }
    let mut union_vec = vec![0u8; alignment_length.div_ceil(8)];
    for gap in gaps {
        bitwise_or_assign(&mut union_vec, gap);
    }
    let gapped_columns = alignment_length - gap_free_columns;
    gapped_columns - count_bits(&union_vec)
}

/// Removes all-gap columns from sequences and filters out excluded sequences.
#[must_use]
pub fn remove_all_gap_columns(
    sequences: &[Vec<u8>],
    headers: &[Vec<u8>],
    excluded: &HashSet<usize>,
) -> (Vec<Vec<u8>>, Vec<Vec<u8>>) {
    if sequences.is_empty() {
        return (Vec::new(), Vec::new());
    }

    let included_indices: Vec<usize> = (0..sequences.len())
        .filter(|idx| !excluded.contains(idx))
        .collect();

    if included_indices.is_empty() {
        return (Vec::new(), Vec::new());
    }

    let seq_len = sequences[0].len();
    let mut gap_columns = vec![true; seq_len];
    for &idx in &included_indices {
        for (pos, &byte) in sequences[idx].iter().enumerate() {
            if !is_gap_char(byte) {
                gap_columns[pos] = false;
            }
        }
    }

    let mut final_sequences = Vec::with_capacity(included_indices.len());
    let mut final_headers = Vec::with_capacity(included_indices.len());

    for &idx in &included_indices {
        let new_seq: Vec<u8> = sequences[idx]
            .iter()
            .enumerate()
            .filter(|&(pos, _)| !gap_columns[pos])
            .map(|(_, &byte)| byte)
            .collect();
        final_sequences.push(new_seq);
        final_headers.push(headers[idx].clone());
    }

    (final_sequences, final_headers)
}
