//! Error types for the `MaxAlign` application.

use std::io;
use std::path::PathBuf;
use thiserror::Error;

/// The primary error type for `MaxAlign` operations.
#[derive(Debug, Error)]
pub enum Error {
    #[error("failed to parse FASTA input: {0}")]
    FastaParse(String),

    #[error("input file is empty")]
    EmptyInput,

    #[error("failed to write output: {0}")]
    WriteOutput(#[from] io::Error),

    #[error("failed to write report to '{path}': {source}")]
    ReportWrite {
        path: PathBuf,
        #[source]
        source: io::Error,
    },

    #[error("failed to write headers list to '{path}': {source}")]
    HeadersListWrite {
        path: PathBuf,
        #[source]
        source: io::Error,
    },
}

pub type Result<T> = std::result::Result<T, Error>;
