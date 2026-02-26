# thecorrelator

A high-performance CLI tool for computing pairwise column correlations from tabular CSV data. For every column pair it reports both **Pearson** and **Spearman** correlation coefficients along with two-tailed p-values.

## Features

- **Two operation modes:**
  - *Within-matrix:* all N*(N-1)/2 column pairs within a single file
  - *Cross-matrix:* every column in file 1 against every column in file 2 (N×M pairs)
- Multi-threaded parallel computation via **OpenMP**
- Reads plain `.csv` and gzip-compressed `.csv.gz` input
- Writes results to stdout or a gzip-compressed `.tsv.gz` file

## Dependencies

- [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/)
- zlib
- A C++11-capable compiler with OpenMP support (e.g. `g++`)

## Build

### CMake (recommended)

```bash
cmake -B build
cmake --build build
# binary is at build/correlate
```

To set the number of parallel compile jobs:

```bash
cmake --build build --parallel 4
```

To install to a custom prefix:

```bash
cmake -B build -DCMAKE_INSTALL_PREFIX=/usr/local
cmake --build build
cmake --install build
```

### Manual (g++)

```bash
g++ -o correlate correlator_parallel.cpp -std=c++11 -lz -fopenmp -lgsl -lgslcblas
```

## Usage

```bash
# Within-matrix: all pairwise correlations among columns in one file
./correlate data.csv [-o results.tsv.gz] [-t <threads>]

# Cross-matrix: correlations between columns of two files
./correlate data1.csv data2.csv [-o results.tsv.gz] [-t <threads>]
```

### Options

| Flag | Description |
|---|---|
| `-o <file>` | Write output to a gzip-compressed TSV file instead of stdout |
| `-t <N>` / `--threads <N>` | Number of OpenMP threads (default: 1) |
| `-h` / `--help` | Print usage information |

## Input Format

Input files must be semicolon-delimited in R's `write.csv2` format:
- First row: column headers
- First column: row names
- Decimal separator: comma (e.g. `3,14` for 3.14)

## Output Format

Tab-separated with columns: `Column1`, `Column2`, `Pearson`, `Pearson_pval`, `Spearman`, `Spearman_pval`.

## License

GNU GPL v3 — see [LICENSE](LICENSE).
