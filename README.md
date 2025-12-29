# maxalign-rs

`maxalign-rs` is a Rust reimplementation of the [MaxAlign](https://services.healthtech.dtu.dk/services/MaxAlign-1.1/) algorithm, offering up to 25-fold speedup relative to the original Perl implementation.

MaxAlign improves multiple sequence alignments by selectively removing sequences that introduce excessive gaps. It maximizes the alignment area (the product of the number of retained sequences and the number of gap-free columns), making more of the alignment usable for downstream analyses such as phylogenetic inference.

## Background

Gaps in sequence alignments can be problematic for some analyses. In phylogenetic inference, for example, columns containing a large number of gaps are frequently discarded entirely. `maxalign-rs` addresses this by identifying sequences that, when removed, yield the largest usable alignment area.

Consider this alignment:

```
            │ 1   2   3   4   5   6   7   8
────────────┼──────────────────────────────
Sequence A  │ -   -   -   A   C   C   T   A
Sequence B  │ -   -   -   A   C   C   T   A
Sequence C  │ G   G   T   A   C   C   G   A
Sequence D  │ G   G   T   A   G   G   -   -
Sequence E  │ C   G   T   A   C   G   G   -
Sequence F  │ T   G   A   A   T   G   C   A
Sequence G  │ T   G   T   A   C   G   C   A
```

Only columns 4, 5, and 6 are gap-free, giving an alignment area of 21 (3 columns ✕ 7 sequences). If sequences A and B are removed, columns 1-6 become gap-free, increasing the alignment area to 30 (6 columns ✕ 5 sequences). MaxAlign automates this optimization.

For more details, refer to the [documentation](https://services.healthtech.dtu.dk/services/MaxAlign-1.1/) of the original implementation or to the associated [publication](https://doi.org/10.1186/1471-2105-8-312).

## Installation

`maxalign-rs` binaries are available for download in the [releases](https://github.com/apcamargo/maxalign-rs/releases/) section of this repository. Alternatively, you can install it from Bioconda using Pixi:

```sh
pixi global install -c bioconda maxalign-rs
```

## Usage

```
maxalign-rs [OPTIONS] [INPUT] [OUTPUT]
```

### Arguments

| Argument | Description | Default |
|----------|-------------|---------|
| `INPUT` | Input FASTA file | `-` (stdin) |
| `OUTPUT` | Output FASTA file | `-` (stdout) |

### Options

| Option | Description | Default |
|--------|-------------|---------|
| `-m`, `--heuristic-method` | Heuristic method: 1 (no synergy), 2 (pairwise synergy), 3 (three-way synergy) | `2` |
| `-i`, `--max-iterations` | Maximum number of iterations (-1 for unlimited) | `-1` |
| `-o`, `--refinement` | Perform refinement using the optimal branch-and-bound algorithm | off |
| `-t`, `--improvement-threshold` | Stop iterating if the relative improvement is below this threshold | `0.0` |
| `-s`, `--excluded-seqs-threshold` | Stop iterating if the fraction of excluded sequences is above this threshold | `1.0` |
| `-k`, `--keep-sequence` | Sequence to always retain (can be specified multiple times) | |
| `-r`, `--report` | Report file path | |
| `--retained-sequences` | Write a list of retained sequences to file | |
| `--excluded-sequences` | Write a list of excluded sequences to file | |
| `-v`, `--verbosity` | Verbosity level (`-v` for normal logging, `-vv` for detailed logging) | off |
| `-h`, `--help` | Print help | |
| `-V`, `--version` | Print version | |

## Examples

### Basic usage

Process an input alignment in the FASTA format and write the optimized result:

```sh
maxalign-rs input.fasta output.fasta
```

Compressed FASTA files in the `.gz`, `.bz2`, `.xz`, and `.zst` formats are supported.

### Using pipes

Read from stdin and write to stdout:

```sh
cat input.fasta | maxalign-rs > output.fasta
```

### Generate a report

Create a detailed Markdown report of the optimization:

```sh
maxalign-rs input.fasta -o output.fasta -r report.md
```

### Limit sequence removal

Stop if more than 20% of sequences would be excluded:

```sh
maxalign-rs input.fasta output.fasta -s 0.2
```

### Protect specific sequences

Keep certain sequences regardless of their gap content. The `-k` option can be specified multiple times to protect multiple sequences:

```sh
maxalign-rs input.fasta output.fasta -k "seq1" -k "seq2" -k "seq3"
```

This ensures that `seq1`, `seq2`, and `seq3` will never be excluded, even if removing them would otherwise improve the alignment area.

### Early stopping

Stop iterating if improvement drops below 1%:

```sh
maxalign-rs input.fasta output.fasta -t 0.01
```

### Export sequence lists

Write lists of included and excluded sequence accessions:

```sh
maxalign-rs input.fasta output.fasta \
    --retained-sequences retained.txt \
    --excluded-sequences excluded.txt
```

### Use a different heuristic method

Use the faster method 1 (no synergy) for very large alignments:

```sh
maxalign-rs input.fasta output.fasta -m 1
```

Or the more thorough method 3 (triple synergy):

```sh
maxalign-rs input.fasta output.fasta -m 3
```

### Ensure optimal solution

Use the branch-and-bound algorithm to refine the results provided by the heuristic algorithm and guarantee an optimal solution:

```sh
maxalign-rs input.fasta output.fasta --refinement
```

Keep in mind that this algorithm performs an exhaustive search and will be very slow for large alignments.

## Citation

If you use `maxalign-rs` in your work, please cite the original paper:

> Gouveia-Oliveira, Rodrigo, Peter W. Sackett, and Anders G. Pedersen. **"MaxAlign: maximizing usable data in an alignment."** *BMC Bioinformatics* 8.1 (2007): 312.
