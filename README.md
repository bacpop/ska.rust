# Split K-mer Analysis (version 2) <img src='ska_logo.png' align="right" height="100" />

<!-- badges: start -->
[![Cargo Build & Test](https://github.com/bacpop/ska.rust/actions/workflows/ci.yml/badge.svg)](https://github.com/bacpop/ska.rust/actions/workflows/ci.yml)
[![docs.rs](https://img.shields.io/docsrs/ska)](https://docs.rs/ska)
[![Clippy check](https://github.com/bacpop/ska.rust/actions/workflows/clippy.yml/badge.svg)](https://github.com/bacpop/ska.rust/actions/workflows/clippy.yml)
[![codecov](https://codecov.io/gh/bacpop/ska.rust/branch/master/graph/badge.svg?token=FZXT39NKA3)](https://codecov.io/gh/bacpop/ska.rust)
[![Crates.io](https://img.shields.io/crates/v/ska)](https://crates.io/crates/ska)
[![GitHub release (latest SemVer)](https://img.shields.io/github/v/release/bacpop/ska.rust)](https://github.com/bacpop/ska.rust/releases)
<!-- badges: end -->

## Description

This is a reimplementation of the [SKA package](https://github.com/simonrharris/SKA)
in the rust language, by Johanna von Wachsmann, Simon Harris and John Lees. We are also grateful to have
received user contributions from:

- Romain Derelle
- Tommi Maklin
- Joel Hellewell
- Timothy Russell
- Nicholas Croucher
- Dan Lu

Split k-mer analysis (version 2) uses exact matching of split k-mer sequences to align closely related sequences, typically small haploid genomes such as bacteria and viruses.

SKA can only align SNPs further than the k-mer length apart, and does not use a gap penalty approach or give alignment scores. But the advantages are speed and flexibility, particularly the ability to run on a reference-free manner (i.e. including accessory genome variation) on both assemblies and reads.

## Documentation

Can be found at https://docs.rs/ska. We also have some tutorials available:

- [From genomes to trees](https://www.bacpop.org/guides/building_trees_with_ska/).

## Installation

Choose from:

1. Download a binary from the [releases](https://github.com/bacpop/ska.rust/releases).
2. Use `cargo install ska` or `cargo add ska`.
3. Use `conda install -c bioconda ska2` (note the two!).
4. Build from source

For 2) or 4) you must have the [rust toolchain](https://www.rust-lang.org/tools/install) installed.

### OS X users

If you have an M1/M2 (arm64) Mac, we aren't currently automatically building binaries, so
would recommend either option 2) or 4) for best performance.

If you get a message saying the binary isn't signed by Apple and can't be run,
use the following command to bypass this:
```
xattr -d "com.apple.quarantine" ./ska
```
### Build from source

1. Clone the repository with `git clone`.
2. Run `cargo install --path .` or `RUSTFLAGS="-C target-cpu=native" cargo install --path .` to optimise for your machine.

## Differences from SKA1

Optimisations include:

- Integer DNA encoding, optimised parsing from FASTA/FASTQ.
- Faster dictionaries.
- Full parallelisation of build phase.
- Smaller, standardised input/output files. Faster to save/load.
- Reduced memory footprint and increased speed with read filtering.

And other improvements:

- IUPAC uncertainty codes for multiple copy split k-mers.
- Uncertainty with self-reverse-complement split k-mers (palindromes).
- Fully dynamic files (merge, delete samples).
- Native VCF output for map.
- Support for known strand sequence (e.g. RNA viruses).
- Stream to STDOUT, or file with `-o`.
- Simpler command line combining `ska fasta`, `ska fastq`, `ska alleles` and `ska merge` into the new `ska build`.
- Option for single commands to run `ska align` or `ska map`.
- New coverage model for filtering FASTQ files with `ska cov`.
- Logging.
- CI testing.

All of which make `ska.rust` run faster and with smaller file size and memory
footprint than the original.

## Planned features

- Sparse data structure which will reduce space and make parallelisation more efficient. [Issue #47](https://github.com/bacpop/ska.rust/issues/47).
- 'fastcall' mode. [Issue #52](https://github.com/bacpop/ska.rust/issues/52).

## Feature ideas (not definitely planned)

- Add support for ambiguity in VCF output (`ska map`). [Issue #5](https://github.com/bacpop/ska.rust/issues/5).
- Non-serial loading of .skf files (for when they are very large). [Issue #22](https://github.com/bacpop/ska.rust/issues/22).
- Alternative mixture models for read error correction. [Issue #50](https://github.com/bacpop/ska.rust/issues/50).

## Things you can no longer do

- Use k > 63 (shouldn't be necessary? Let us know if you need this and why).
- `ska annotate` (use bedtools).
- `ska compare`, `ska humanise`, `ska info` or `ska summary` (replaced by `ska nk --full-info`).
- `ska unique` (you can parse `ska nk --full-info` if you want this functionality, but we didn't think it's used much).
- `ska type` (use [PopPUNK](https://github.com/bacpop/PopPUNK) instead of MLST 🙂)
- Ns are always skipped, and will not be found in any split k-mers.
- `.skf` files are not backwards compatible with version 1.
