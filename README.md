# Split K-mer Analysis (version 2)

<!-- badges: start -->
[![Cargo Build & Test](https://github.com/bacpop/ska.rust/actions/workflows/ci.yml/badge.svg)](https://github.com/bacpop/ska.rust/actions/workflows/ci.yml)<!-- badges: end -->

## Installation

Choose from:

1. Download a binary from the [releases](https://github.com/bacpop/ska.rust/releases).
2. Use `cargo install ska`.
3. Use `conda install -c bioconda ska`.
4. Build from source

For 2) or 4) you must have the [rust toolchain](https://www.rust-lang.org/tools/install) installed.
### Build from source

1. Clone the repository with `git clone`.
2. Run `cargo install --path .` or `RUSTFLAGS="-C target-cpu=native" cargo install --path .` to optimise for your machine.

## Documentation

Can be found at https://docs.rs/ska.

## Description

This is a reimplementation of Simon Harris' [SKA package](https://github.com/simonrharris/SKA)
in the rust language, by Johanna von Wachsmann, Simon Harris and John Lees.

> SKA (Split Kmer Analysis) is a toolkit for prokaryotic (and any other small, haploid) DNA sequence analysis using split kmers. A split kmer is a pair of kmers in a DNA sequence that are separated by a single base. Split kmers allow rapid comparison and alignment of small genomes, and is particulalry suited for surveillance or outbreak investigation. SKA can produce split kmer files from fasta format assemblies or directly from fastq format read sequences, cluster them, align them with or without a reference sequence and provide various comparison and summary statistics. Currently all testing has been carried out on high-quality Illumina read data, so results for other platforms may vary.

Optimisations include:

- Integer DNA encoding, optimised parsing from FASTA/FASTQ.
- Faster dictionaries.
- Full parallelisation of build phase.
- Smaller, standardised input/output files.
- Reduced memory footprint with read filtering.

And other improvements:

- IUPAC uncertainty codes for multiple copy k-mers.
- Fully dynamic files (merge, delete samples).
- Native VCF output for map.
- Support for known strand sequence (e.g. RNA viruses).
- Stream to STDOUT, or file with `-o`.
- Simpler command line combining `ska fasta`, `ska fastq`, `ska alleles` and `ska merge` into the new `ska build`.
- Option for single commands to run `ska align` or `ska map`.
- Logging.
- CI testing.

All of which make ska.rust run faster and with smaller file size and memory
footprint than the original.

Things you can no longer do:

- Use k > 31 (shouldn't be necessary).
- `ska annotate` (use bedtools).
- `ska compare`, `ska humanise`, `ska info` or `ska summary` (use `ska nk --full-info`).
- `ska distance` and `ska unique` (use [pp-sketchlib](https://github.com/bacpop/pp-sketchlib)).
- `ska type` (use [PopPUNK](https://github.com/bacpop/PopPUNK) instead of MLST 🙂)
- Ns are always skipped, and will not be found in any split k-mers.
- `.skf` files are not backwards compatible with version 1.