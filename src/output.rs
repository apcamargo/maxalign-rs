//! Output utilities for writing FASTA files and header lists.

use crate::error::{Error, Result};
use crate::fasta::get_record_accession_string;
use clio::Output;
use itertools::Itertools;
use std::collections::HashSet;
use std::io::{BufWriter, Write};
use std::path::Path;

const FASTA_LINE_WIDTH: usize = 80;

/// Writes sequences in FASTA format to the given output.
pub fn write_fasta(sequences: &[Vec<u8>], headers: &[Vec<u8>], output: &mut Output) -> Result<()> {
    if sequences.is_empty() {
        return Ok(());
    }

    for (header, seq) in headers.iter().zip_eq(sequences) {
        writeln!(output, ">{}", String::from_utf8_lossy(header))?;
        for chunk in seq.chunks(FASTA_LINE_WIDTH) {
            writeln!(output, "{}", String::from_utf8_lossy(chunk))?;
        }
    }

    Ok(())
}

/// Writes a list of headers to a file (included or excluded based on the flag).
pub fn write_headers_list(
    path: impl AsRef<Path>,
    headers: &[Vec<u8>],
    excluded: &HashSet<usize>,
    included: bool,
) -> Result<()> {
    let path = path.as_ref();
    let file = std::fs::File::create(path).map_err(|e| Error::HeadersListWrite {
        path: path.to_path_buf(),
        source: e,
    })?;
    let mut writer = BufWriter::new(file);

    for (idx, header) in headers.iter().enumerate() {
        if excluded.contains(&idx) != included {
            let accession = get_record_accession_string(header).unwrap_or_default();
            writeln!(writer, "{}", accession).map_err(|e| Error::HeadersListWrite {
                path: path.to_path_buf(),
                source: e,
            })?;
        }
    }

    writer.flush().map_err(|e| Error::HeadersListWrite {
        path: path.to_path_buf(),
        source: e,
    })?;

    Ok(())
}
