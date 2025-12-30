# maxalign-rs

`maxalign-rs` is a Rust reimplementation of the [MaxAlign](https://services.healthtech.dtu.dk/services/MaxAlign-1.1/) algorithm, offering up to 25-fold speedup relative to the original Perl implementation.

MaxAlign improves multiple sequence alignments by selectively removing sequences that introduce excessive gaps. It maximizes the alignment area (the product of the number of retained sequences and the number of gap-free columns), making more of the alignment usable for downstream analyses such as phylogenetic inference.

## Background

Gaps in sequence alignments can be problematic for some analyses. In phylogenetic inference, for example, columns containing a large number of gaps are frequently discarded entirely. `maxalign-rs` addresses this by identifying sequences that, when removed, yield the largest usable alignment area.

Consider this alignment:

```
            │ 1   2   3   4   5   6   7   8
────────────┼──────────────────────────────
Sequence A  │ -   -   -   A   C   C   t   a
Sequence B  │ -   -   -   A   C   C   t   a
Sequence C  │ g   g   t   A   C   C   g   a
Sequence D  │ g   g   t   A   G   G   -   -
Sequence E  │ c   g   t   A   C   G   g   -
Sequence F  │ t   g   a   A   T   G   c   a
Sequence G  │ t   g   t   A   C   G   c   a
```

Only columns 4, 5, and 6 are gap-free, giving an alignment area of 21 (3 columns × 7 sequences). If sequences A and B are removed, columns 1-6 become gap-free, increasing the alignment area to 30 (6 columns × 5 sequences).

```
            │ 1   2   3   4   5   6   7   8
────────────┼──────────────────────────────
Sequence C  │ G   G   T   A   C   C   g   a
Sequence D  │ G   G   T   A   G   G   -   -
Sequence E  │ C   G   T   A   C   G   g   -
Sequence F  │ T   G   A   A   T   G   c   a
Sequence G  │ T   G   T   A   C   G   c   a
```

MaxAlign automates this optimization process by iteratively identifying and removing sequences whose exclusion maximizes the resulting alignment area. For more details, refer to the [documentation](https://services.healthtech.dtu.dk/services/MaxAlign-1.1/) of the original implementation or to the associated [publication](https://doi.org/10.1186/1471-2105-8-312).

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
| `-o`, `--refinement` | Perform refinement using the branch-and-bound algorithm to find the optimal solution | off |
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
# Reading from a file and writing to a file
maxalign-rs input.fasta output.fasta
# Reading from stdin and writing to stdout
cat input.fasta | maxalign-rs > output.fasta
```

`maxalign-rs` accepts compressed input FASTA files in the `.gz`, `.bz2`, `.xz`, and `.zst` formats.

### Use a different heuristic method

MaxAlign applies a greedy heuristic that iteratively removes sequences to maximize the alignment area, defined as the number of retained sequences multiplied by the number of gap-free columns. Candidate removals are derived from gap patterns, and the process stops when no further improvement is possible.

Three heuristic methods are available, differing in how they account for interactions between sequence removals:

- **Method 1 (no synergy):** Evaluates each candidate sequence set removal independently. For example, if removing sequences ${A, B}$ would free 3 columns and removing ${D, E}$ would free 2 columns, it considers these independently. Fastest option, but may miss improvements that require removing multiple sequences together. Recommended for very large alignments.
- **Method 2 (pairwise synergy, default):** Considers that removing two sequence sets together might free more columns than the sum of removing them separately. For instance, removing ${A, B}$ frees 3 columns and ${D, E}$ frees 2 columns, but removing both together might free more than 5 columns. Balances runtime and solution quality and is recommended for most use cases.
- **Method 3 (three-way synergy):** Extends the logic to three sequence sets, checking whether removing three sets together provides additional benefit beyond pairwise combinations. Can yield marginally better results at the cost of increased computation time.

In `maxalign-rs`, you can select the heuristic method using the `-m` option:

```sh
# Use the faster method 1 (no synergy) for very large alignments
maxalign-rs input.fasta output.fasta -m 1
# Use method 3 (three-way synergy) for potentially better results (higher alignment area)
maxalign-rs input.fasta output.fasta -m 3
```

### Ensure maximal alignment area

While the heuristic algorithms are fast and find the optimal solution in most cases, you can use the branch-and-bound algorithm to guarantee finding the absolute best solution. The branch-and-bound algorithm exhaustively searches through all possible combinations of sequence removals to find the truly optimal solution.

When `-o` is specified, `maxalign-rs` first runs the heuristic algorithm to get a good solution quickly, then applies the branch-and-bound algorithm to verify whether this solution is optimal or to find a better one if it exists:

```sh
maxalign-rs input.fasta output.fasta -o
```

Keep in mind that this algorithm performs an exhaustive search and will be very slow for large alignments.

### Limit sequence removal

You can limit the number of sequences removed during optimization by stopping the process early based on two criteria: the fraction of sequences excluded from the alignment (`-s`), and the relative improvement in alignment area between iterations (`-t`). The process stops early if the excluded fraction exceeds a specified threshold or if the relative improvement falls below a specified threshold.

```sh
# Stops iterating if more than 20% of sequences would be excluded
maxalign-rs input.fasta output.fasta -s 0.2
# Stop iterating if relative improvement between iterations is less than 1%
maxalign-rs input.fasta output.fasta -t 0.01
```

### Protect sequences from removal

`maxalign-rs` allows you to force specific sequences to be retained in the output alignment, even if their removal would increase the alignment area. Use the `-k` option to specify sequences to protect. This option can be provided multiple times:

```sh
maxalign-rs input.fasta output.fasta -k "seq1" -k "seq2" -k "seq3"
```

This ensures that `seq1`, `seq2`, and `seq3` will always be included in the final alignment.

### Generate a report

Generate a detailed Markdown report summarizing the optimization process, including the number of retained and excluded sequences, changes in alignment area across iterations, and the final optimization outcome:

```sh
maxalign-rs input.fasta output.fasta -r report.md
```

### Export sequence lists

Write plain-text files that record which sequence accessions were retained and which were excluded during the alignment optimization procedure:

```sh
maxalign-rs input.fasta output.fasta \
    --retained-sequences retained.txt \
    --excluded-sequences excluded.txt
```

## Citation

If you use `maxalign-rs` in your work, please cite the original paper:

> Gouveia-Oliveira, Rodrigo, Peter W. Sackett, and Anders G. Pedersen. [**"MaxAlign: maximizing usable data in an alignment"**](https://doi.org/10.1186/1471-2105-8-312). *BMC Bioinformatics* 8.1 (2007): 312.
