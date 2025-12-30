//! Bit manipulation utilities for efficient set operations.

const BITS_PER_BYTE: usize = 8;

/// Counts the number of set bits (1s) in a byte slice.
#[must_use]
pub fn count_bits(bytes: &[u8]) -> usize {
    bytes.iter().map(|&b| b.count_ones() as usize).sum()
}

/// Computes the bitwise OR (union) of two byte slices.
#[must_use]
pub fn bitwise_or(a: &[u8], b: &[u8]) -> Vec<u8> {
    a.iter().zip(b).map(|(&x, &y)| x | y).collect()
}

/// Computes the bitwise OR of two byte slices in place.
pub fn bitwise_or_assign(dest: &mut [u8], src: &[u8]) {
    for (d, &s) in dest.iter_mut().zip(src) {
        *d |= s;
    }
}

/// Sets a single bit at the specified position.
pub fn set_bit(vec: &mut [u8], position: usize) {
    let byte_index = position / BITS_PER_BYTE;
    let bit_offset = position % BITS_PER_BYTE;
    if byte_index < vec.len() {
        vec[byte_index] |= 1 << bit_offset;
    }
}

/// Packs a slice of booleans into a bit-packed byte vector.
#[must_use]
pub fn pack_bools_to_bits(bools: &[bool]) -> Vec<u8> {
    bools
        .chunks(BITS_PER_BYTE)
        .map(|chunk| {
            let mut byte = 0u8;
            for (i, &b) in chunk.iter().enumerate() {
                if b {
                    byte |= 1 << i;
                }
            }
            byte
        })
        .collect()
}

/// Returns the indices of all set bits in a bit-packed byte vector.
#[must_use]
pub fn get_set_bit_indices(bytes: &[u8], count: usize) -> Vec<usize> {
    let mut indices = Vec::with_capacity(count_bits(bytes));
    for (byte_idx, &byte) in bytes.iter().enumerate() {
        if byte == 0 {
            continue;
        }
        for bit_idx in 0..BITS_PER_BYTE {
            let i = byte_idx * BITS_PER_BYTE + bit_idx;
            if i >= count {
                break;
            }
            if (byte >> bit_idx) & 1 == 1 {
                indices.push(i);
            }
        }
    }
    indices
}

/// Computes the population count of the bitwise OR of two byte slices.
#[must_use]
pub fn count_bits_union(a: &[u8], b: &[u8]) -> usize {
    a.iter()
        .zip(b)
        .map(|(&x, &y)| (x | y).count_ones() as usize)
        .sum()
}

/// Computes the population count of the bitwise OR of three byte slices.
#[must_use]
pub fn count_bits_union_triple(a: &[u8], b: &[u8], c: &[u8]) -> usize {
    a.iter()
        .zip(b)
        .zip(c)
        .map(|((&x, &y), &z)| (x | y | z).count_ones() as usize)
        .sum()
}
