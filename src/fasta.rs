//! FASTA file parsing utilities.

use crate::error::{Error, Result};
use clio::Input;
use itertools::Itertools;
use log::warn;
use needletail::{parse_fastx_file, parse_fastx_stdin};
use std::collections::HashSet;

/// Extracts the accession (first word) from a FASTA header.
pub fn get_record_accession_string(record_header: &[u8]) -> Option<String> {
    let accession = record_header
        .split(|&b| matches!(b, b' ' | b'\t' | b'\n' | b'\x0C' | b'\r'))
        .next();
    match accession {
        Some(acc) if !acc.is_empty() => Some(String::from_utf8_lossy(acc).into_owned()),
        _ => None,
    }
}

/// Parsed sequence data from a FASTA file.
pub struct SequenceData {
    pub headers: Vec<Vec<u8>>,
    pub sequences: Vec<Vec<u8>>,
    pub longest_length: usize,
    pub keep_indices: HashSet<usize>,
}

/// Parses a FASTA file and returns the sequence data.
pub fn parse_fasta(input: &Input, keep_sequence: &[String]) -> Result<SequenceData> {
    let reader = if input.is_std() {
        parse_fastx_stdin()
    } else {
        if input.is_empty().unwrap_or(false) {
            return Err(Error::EmptyInput);
        }
        parse_fastx_file(input.path().to_path_buf())
    };

    let mut reader = reader.map_err(|e| Error::FastaParse(e.to_string()))?;

    let mut headers = Vec::new();
    let mut sequences = Vec::new();
    let mut longest_length = 0;
    let keep_set: HashSet<&str> = keep_sequence.iter().map(String::as_str).collect();

    while let Some(record) = reader.next() {
        let record = record.map_err(|e| Error::FastaParse(e.to_string()))?;

        let header_bytes = record.id().to_vec();
        let mut sequence_bytes = record.seq().to_vec();
        sequence_bytes.retain(|&b| !b.is_ascii_whitespace());

        longest_length = longest_length.max(sequence_bytes.len());

        headers.push(header_bytes);
        sequences.push(sequence_bytes);
    }

    if sequences.is_empty() {
        return Err(Error::EmptyInput);
    }

    let (min_length, longest_length_found) = sequences
        .iter()
        .map(Vec::len)
        .minmax()
        .into_option()
        .unwrap_or((0, 0));

    if min_length != longest_length_found {
        warn!(
            "Sequences have different lengths ({min_length} to {longest_length_found}). Shorter sequences will be padded with gaps."
        );
    }

    let mut keep_indices = HashSet::new();
    let mut found_keep_sequence: HashSet<String> = HashSet::new();

    for (idx, header) in headers.iter().enumerate() {
        if let Some(accession) = get_record_accession_string(header)
            && keep_set.contains(accession.as_str())
        {
            keep_indices.insert(idx);
            found_keep_sequence.insert(accession);
        }
    }

    for accession in keep_sequence {
        if !found_keep_sequence.contains(accession) {
            warn!(
                "Must-retain sequence '{}' was not found in the input alignment",
                accession
            );
        }
    }

    Ok(SequenceData {
        headers,
        sequences,
        longest_length,
        keep_indices,
    })
}
