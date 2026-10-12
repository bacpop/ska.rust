//! Split k-mer analysis (version 2) uses exact matching of split k-mer sequences to align closely related
//! sequences, typically small haploid genomes such as bacteria and viruses.
//!
//! SKA can only align SNPs further than the k-mer length apart,
//! and does not use a gap penalty approach or give alignment scores.
//! But the advantages are speed and flexibility, particularly the ability to
//! run on a reference-free manner (i.e. including accessory genome variation)
//! on both assemblies and reads.
//!
//! ## Details
//!
//! A split k-mer is a sequence of bases with the middle position removed. For
//! example, the split 9-mer pattern is `XXXX-XXXX`, and the split 9-mers of the
//! sequence `ACGAGAGTCTT` are:
//!
//! | Split k-mer | Middle base |
//! |-------------|-------------|
//! | `ACGAAGTC`  | `G`         |
//! | `CGAGGTCT`  | `A`         |
//! | `GAGATCTT`  | `G`         |
//!
//! Which is broadly how SKA represents sequences in `.skf` files, as a dictionary with split
//! k-mers as keys and the middle base as values. These dictionaries are merged
//! by exactly matching the split k-mers, such that the middle positions are aligned (but
//! unordered). Split k-mers can also be matched to those from a reference sequence
//! to give an ordered pseudoalignment.
//!
//! Various optimisations are used to make this as fast as possible. For a more thorough comparison with version 1.0 of SKA, see the
//! [github description](https://github.com/bacpop/ska.rust/blob/master/README.md).
//!
//! *NB As split k-mers are even lengths, it is possible that they are their
//! own reverse complement. The original version of ska would randomly pick a strand,
//! possibly introducing a SNP across samples. This version uses an ambiguous middle
//! base (W for A/T; S for C/G) to represent this case.*
//!
//! Tutorials:
//! - [From genomes to trees](https://www.bacpop.org/guides/building_trees_with_ska/).
//! - [Filtering options](https://www.bacpop.org/guides/snp_alignment_with_ska/).
//!
//! Papers:
//! - [Seamless, rapid, and accurate analyses of outbreak genomic data using split k-mer analysis ](https://genome.cshlp.org/content/34/10/1661.abstract).
//! - [skalo: using SKA split k-mers with coloured de Brujin graphs to genotype indels](https://academic.oup.com/mbe/article/42/4/msaf077/8103706).
//!
//! Command line usage follows. For API documentation and usage, see the [end of this section](#api-usage).
//!
//! # Usage
//!
//! `.skf` files represent merged split k-mers from multiple sequence files. They
//! are created with `ska build`. You can subsequently use `ska align` to write
//! out an unordered alignment from these files, or `ska map` to write an alignment
//! ordered versus a reference sequence.
//!
//! Alternatively, both `ska align` and `ska map` can take input sequence directly to obtain output alignments
//! in a single command and without saving an `.skf` file. NB: This uses the default options
//! of `ska build`, so to change these you will need to run the alignment in two steps.
//!
//! Output from `ska align` and `ska map` is to STDOUT, so you can use a redirect `>` to save to a file or pipe `|`
//! to stream into another compatible program on STDIN. You can also add an output
//! file prefix directly with `-o` (for `ska build` this is required).
//!
//! ## Important notes
//!
//! - In this version of ska the k-mer length is the _total_ split k-mer size,
//!   whereas in ska1 the k-mer length is half the split length. So the default of
//!   k=15 in ska1 corresponds to choosing `-k 31` in this version.
//! - If you are using FASTQ input files (reads), these must be provided as two
//!   deinterleaved files using the `-f` option to `ska build`. Providing them as
//!   a single file will treat them as FASTA, and not correct sequencing errors.
//! - If you are running `ska map`, it is more memory efficient to run these one at
//!   a time, rather than merging a single `.skf` file. A command for doing this
//!   in parallel is listed below.
//!
//! ## Common options
//!
//! Version can be viewed by running `ska -V`.
//!
//! Details and progress messages are written on STDERR. You can see more logging
//! information by adding the verbose flag `-v`.
//!
//! ## ska build
//!
//! This command creates an `.skf` file from sequence (.fasta/.fasta.gz/.fastq/.fastq.gz) input.
//! K-mer size must be odd, greater than 5, and a maximum of 63. Smaller k-mers
//! are more sensitive and can find closer positions, but are less specific so
//! may lead to more repeated split k-mers and ambiguous bases. Using k <= 31 uses
//! 64-bit integers and may be faster than 31 < k <= 63, which uses 128-bits.
//!
//! This is typically the most computationally intensive step of `ska`, and
//! multiple `--threads` can be used to split the work over multiple CPU cores.
//!
//! Build from two input FASTA files with a k-mer size of 31:
//! ```bash
//! ska build -o seqs -k 31 assemblies/seq1.fa assemblies/seq2.fa
//! ```
//! This will assume sample names of `seq1` and `seq2` –- the base path and with known fastx
//! file extensions removed. If you know the strand,
//! for example when exclusively using reference sequences or single stranded viruses, add `--single-strand`
//! to ignore reverse complements.
//!
//! ### Read files
//!
//! To use FASTQ files, or specify sample names, or more easily input a larger number of input files,
//! you can provide a tab separated file list via `-f` instead of listing files. For example
//! a file called `input_sequence.txt` to describe `sample_1` (an assembly) and `sample_2` (paired reads):
//! ```text
//! sample_1    assemblies/seq1.fa
//! sample_2    reads/seq2_1.fastq.gz   reads/seq2_2.fastq.gz
//! ```
//! You'd run:
//! ```bash
//! ska build -o seqs -f input_sequence.txt
//! ```
//! Options for filtering/error correction are:
//! - `--min-count`. Specify a minimum number of appearances a k-mer must have
//!   to be included. This is an effective way of filtering sequencing errors if set
//!   to at least three, but higher may be appropriate for higher coverage data.
//!   A two-step blocked bloom and hash table filter is used for memory efficiency.
//! - `--proportion-reads`. Specify a proportion of reads to use to build the .skf file.
//!   The value of this parameter must be between 0 and 1. If you have high coverage samples
//!   this can be used to reduce the build time.
//! - `--qual-filter`. `none` do not filter based on quality scores.
//!   `middle` (default) filter k-mers where the middle base is below the minimum quality.
//!   `strict` filter k-mers where any base is below the minimum quality.
//! - `--min-qual`. Specify a minimum PHRED score to use in the filter.
//!
//! FASTQ files must be paired end. If you'd like to request more flexibility in
//! this regard please contact us.
//!
//! ### Threads
//!
//! The maximum threads actually used will be a power of two, so if you provided
//! `--threads 6` only four would be used. Additionally, at least ten samples per
//! thread are required so maximums are:
//!
//! | Samples | Maximum threads |
//! |---------|-----------------|
//! | 1-19    | 1               |
//! | 20-39   | 2               |
//! | 40-79   | 4               |
//!
//! and so on. Use `-v` to see a message with the number being used.
//!
//! Using more threads will increase memory use.
//!
//! You can also run blocks of samples independently (e.g. with snakemake or a
//! job array) then use `ska merge` to combine results.
//!
//! ## ska align
//!
//! Creates an alignment from a `.skf` file or sequence files. Sites (columns) are
//! in an arbitrary order. Two basic filters are available: `--min-freq` which
//! sets the maximum number of missing sites; `--filter` which can be set to
//! one of three options:
//! * `no-filter` -- no extra filtering
//! * `no-const` -- no constant sites
//! * `no-ambig-or-const` -- no constant sites, or sites where the only variable base is ambiguous
//!
//! `no-filter` may be useful for ascertainment bias correction in phylogenetic algorithms,
//! but note the flanking split k-mers will never be included. `no-const` is the default.
//! `no-ambig-or-const` can be used when you want to treat any ambiguity as an `N`.
//!
//! With an `.skf` file from `ska build`, constant sites, and no missing variants:
//! ```bash
//! ska align --min-freq 1 --filter no-filter -o seqs seqs.skf
//! ```
//!
//! Another example: directly from FASTA files, with default `build` and `align` settings,
//! and gzipping the output alignment
//! ```bash
//! ska align seq*.fa --threads 8 | gzip -c - > seqs.aln.gz
//! ```
//!
//! ## ska map
//!
//! Creates an alignment from a `.skf` file or sequence files by mapping its
//! split k-mers to split k-mers of a linear reference sequence. This produces pseudoalignment
//! with the same sites/columns as the reference sequence. Sites not mapped will
//! be output as missing '-'.
//!
//! Pass the FASTA file of the reference as the first argument, and the skf file as
//! the second argument:
//! ```bash
//! ska map ref.fa seqs.skf -o ref_mapped.aln
//! ```
//! Add `--repeat-mask` to mask any repeated split k-mers in the reference with 'N'.
//!
//! You can also get a VCF as output, which has rows as variants, and only has the
//! variable sites (but will include unmapped bases as missing).
//! An example command, also demonstrating that everything can be done from input sequences in a single command:
//! ```bash
//! ska map ref.fa seq1.fa seq2.fa -f vcf --threads 2 | bgzip -c - > seqs.vcf.gz
//! ```
//!
//! ## ska lo
//!
//! Converts split k-mers from a `.skf` file into a colored De Bruijn graph and infers indels from graph bubbles and SNPs from variant groups in
//! reference-free mode (as with `ska align`). SNPs are only composed of ATGC variants (no ambigous nucleotides). The same
//! filtering applies to indels.
//! Multithreading ('-t' argument) is not fully optimised - it usually takes 4 threads to halve runtimes.
//!
//! To generate a SNP alignment and an indel VCF file (here named 'test_snps.fas' and 'test_indels.vcf'):
//! ```bash
//! ska lo seqs.skf test
//! ```
//!
//! It can also position SNPs on a reference genome if provided using the '-r' argument. The reference genome should
//! be in FASTA format and composed of a unique sequence. skalo lo will then generate, in addition to the SNP alignment and the indel
//! VCF file, a SNP VCF file and a pseudo-genome alignment (as with `ska map`) that can be used for recombination analyses.
//! In such use case, we recommmend to increase the maximum proportion of missing data allowed per variant ('-m' argument), but not above 0.5.
//! ```bash
//! ska lo seqs.skf test -r reference.fas -m 0.4
//! ```
//! Please note that at the moment indels cannot be positioned on a reference genome.
//!
//!
//! ### Efficiency
//!
//! As `ska map` is independent for each input file, the most efficient way to
//! run this on a number of samples is with gnu parallel, for example:
//! ```bash
//! cat ${FILE_LIST} | parallel -j 8 "ID=\$(basename {} | sed 's/_a.*fasta/_temp_skf/g');
//! ska build -o \$ID -k 31 {} ;
//! MAPOUT=\$(basename {} | sed 's/_a.*fasta/_aln/g');
//!  ska map $REFERENCE_LOC \${ID}.skf -o \${MAPOUT}.aln"
//! cat *.aln > ska_map.aln
//! ```
//!
//! ## ska distance
//!
//! Use to calculate distances between all samples within an `.skf` file. The output
//! will contain the number of SNP differences between all pairs of samples, as
//! well as the proportion of matching k-mers.
//!
//! ```bash
//! ska distance -o distances.txt seqs.skf
//! ```
//!
//! Mismatches is an estimate of the alignable genome fraction based on the unfiltered
//! k-mers.
//!
//! Also consider ambiguous bases by adding `--allow-ambiguous` flag, which will give partial SNP distances.
//! Note that ambiguous bases may overestimate distances due to repeat k-mers. Use `--min-freq` to
//! ignore k-mers only found in some samples at the population level (default = 0.0). F
//! or finer control over filtering, first run `ska weed`
//! on the input .skf.
//!
//! Multiple threads can be used, but this will only be effective with large numbers of samples.
//!
//! The companion script in `scripts/cluster_dists.py` (requires `networkx`) can
//! be used to make single linkage clusters from these distances at given thresholds,
//! and create a visualisation in [Microreact](https://microreact.org/):
//! ```bash
//! ska distance seqs.skf > distances.txt
//! python scripts/cluster_dists.py distances.txt --snps 20 --mismatches 0.05
//! ```
//! If you install `rapidnj` and `biopython` you can also draw an NJ tree from
//! these distances, which will be displayed in Microreact. Use your `--api-key`
//! to directly upload and get the URL printed on the terminal.
//!
//! ## ska merge
//!
//! Use to combine multiple `.skf` files into one, typically for subsequent use in `align` or `map`.
//! This may be particularly useful if `ska build` was run on multiple input files
//! one at a time (for example in a job array).
//!
//! ```bash
//! ska merge -o all_samples sample1.skf sample2.skf sample3.skf
//! ```
//!
//! ## ska delete
//!
//! Remove samples by name from an `.skf` file. All sample names must exist in the
//! file, or an error will be thrown. After removing samples, the input `.skf` file will be overwritten,
//! unless you provide a new one with `-o`.
//! Any split k-mers which have no observed missing bases after removal will also be deleted.
//!
//! ```bash
//! ska delete --skf-file all_samples.skf sample_1 sample_3
//! ```
//! If you wish to delete many samples, you can use `-f` to provide their names
//! by line via an input file.
//!
//! ## ska weed
//!
//! Remove (weed) split k-mers from an `.skf` file. The split k-mers to be removed
//! are generated from a FASTA file, which may for example contain known transposons or
//! other repeat sequences. As with `delete`, the input `.skf` file is overwritten by default.
//! Use `-o` to write the results to a new file instead of overwriting the input `.skf`.
//!
//! ```bash
//! ska weed all_samples.skf MGEs.fa
//! ```
//!
//! In addition, you can also use `ska weed` to filter split k-mers by proportion of samples
//! they appear in, and constant or amibguous middle cases, which does
//! not require a list of split k-mers to be removed (but both can be used together, if you wish).
//!
//! For example, you may want to reduce the size of an `.skf` file for online use.
//! You can do this by removing any constant sites or ambiguous-only s0tes (which are typically unused), and by hard-filtering
//! by frequency (i.e. before writing output):
//! ```bash
//! ska weed --filter no-ambig-or-const --min-freq 0.9 all_samples.skf
//! ```
//!
//! ## ska nk
//! Return information on the *n*umber of *k*-mers in an `.skf` file. This will
//! print on STDOUT:
//! - The version of ska used to build the file.
//! - The k-mer size.
//! - Number of bits used for the split k-mer (64 or 128).
//! - Whether reverse complements were used.
//! - The total number of split k-mers.
//! - The total number of samples.
//! - The sample names.
//! - The number of split k-mers found in each sample.
//!
//! ```bash
//! ska nk all_samples.skf
//! ska_version=0.3.1
//! k=21
//! k_bits=64
//! rc=true
//! k-mers=3228084
//! samples=28
//! sample_names=["19183_4#73", "12673_8#29", "12673_8#31", "12754_5#61", "12754_5#89", "12754_5#47", "12754_5#32", "12754_5#78", "12754_4#85", "12754_5#74", "19183_4#57", "12754_5#36", "19183_4#49", "19183_4#79", "19183_4#60", "12754_5#24", "12754_5#22", "12754_5#71", "12673_8#26", "12754_5#95", "12754_5#86", "12673_8#24", "19183_4#61", "12673_8#41", "12754_4#66", "12754_5#80", "12754_5#84", "19183_4#78"]
//! sample_kmers=[2872587, 2997448, 2949719, 2997496, 2997178, 2912749, 2996491, 2997221, 2949102, 2997454, 2914109, 2912237, 2872518, 2869957, 2872470, 2997992, 2997647, 2958512, 2998099, 2997290, 2950253, 3027707, 2997881, 2907920, 2911447, 2997644, 2944830, 2915080]
//! ```
//!
//! If you add `--full-info`, the split k-mer dictionary will be decoded and printed
//! to STDOUT after the baseline information, for example:
//! ```bash
//! TAAATATC        TAAACGCC        A,-
//! AGACTCTC        TACAGCTA        G,G
//! AAACCATC        AAACACTC        A,-
//! TTAAAAGA        TCTCGTAC        C,C
//! ```
//! This is tab-separated, showing the first and second half of the split k-mer
//! and the middle bases of each sample (comma seperated).
//!
//! NB: Only one split k-mer is shown even if the reverse complement was used.
//! These are not precisely canonical k-mers, as the encoding order `{A, C, T, G}` is used internally.
//! But if you can't find a sequence in your input file, you will find its reverse complement.
//!
//! ## ska cov
//!
//! Estimate a coverage cutoff for use with read data. This will count the split
//! k-mers in a pair of FASTQ samples, and create a histogram of these counts.
//! A mixture model is then fitted to this histogram using maximum likelihood,
//! which can give a suggested cutoff with noisy data.
//!
//! The cutoff will be printed to STDERR. Use `-v` to get additional information on the
//! optimisation progress and result. A table of counts and the fit will be printed
//! to STDOUT, which can then be plotted by the companion script in
//! `scripts/plot_cov.py` (requires `matplotlib`):
//! ```bash
//! ska cov reads_1.fastq.gz reads_2.fastq.gz -k 31 -v > cov_plot.txt
//! python scripts/plot_cov.py cov_plot.txt
//! ```
//!
//! The cutoff can be used with the `--min-count` parameter of `ska build`. For
//! a set of sequence experiments with similar characteristics it may be sufficient
//! to use the same cutoff. Alternatively `ska cov` can be run on every sample
//! independently (`gnu parallel` would be an efficient way to do this).
//!
//! ## ska lo
//!
//! Runs "skalo", a graph-based algorithm for inferring variants from WGS outbreak data.
//! Skalo takes as input a `.skf` file generated by the `ska build` command.
//!
//! To run skalo on a file named `my_file.skf` use the following command `ska lo my_file.skf skalo_results`.
//!
//! The skalo method
//! can be modified using the following flags:
//! * `-r` the path to a reference genome that the detected SNPs can be mapped onto
//! * `-m` the maximum fraction of missing data permissible for a site where a possible SNP has been detected (default = 0.2)
//! * `-d` the depth of the recursive path allowed during graph traversal (default = 4)
//! * `-n` maximum number of internal indel k-mers (default = 2)
//! * `-t` number of threads to use for parallel computation (default = 1)
//!
//!
//! # API usage
//!
//! See the submodule documentation linked below.
//!
//! To use `ska build` in other rust code:
//! ```rust
//! use ska::merge_ska_dict::{InputFastx, build_and_merge};
//! use ska::merge_ska_array::MergeSkaArray;
//! use ska::{QualOpts, QualFilter};
//!
//! // Build, merge
//! let input_files: [InputFastx; 2] = [("test1".to_string(),
//!                                      "tests/test_files_in/test_1.fa".to_string(),
//!                                      None),
//!                                     ("test2".to_string(),
//!                                      "tests/test_files_in/test_2.fa".to_string(),
//!                                      None)];
//! let rc = true;
//! let k = 31;
//! let quality = QualOpts {min_count: 1, min_qual: 0, qual_filter: QualFilter::NoFilter};
//! let threads = 2;
//! let proportion_reads = None;
//! // NB u64 for k<=31, u128 for k<=63
//! let merged_dict =
//!     build_and_merge::<u64>(&input_files, k, rc, &quality, threads, proportion_reads);
//!
//! // Save
//! let ska_array = MergeSkaArray::new(&merged_dict);
//! ska_array
//!     .save(format!("merged.skf").as_str())
//!     .expect("Failed to save output file");
//! ```
//!
//! To use `ska align` in other rust code:
//! ```rust
//! use ska::io_utils::{load_array, set_ostream};
//! use ska::cli::FilterType;
//!
//! // Load an .skf file
//! let threads = 4;
//! let input = vec!["tests/test_files_in/merge.skf".to_string()];
//! let mut ska_array = load_array::<u64>(&input, threads).expect("Could not open input as u64");
//!
//! // Apply filters
//! let min_count = 2;
//! let filter_ambig_as_missing = false;
//! let update_kmers = false;
//! let filter = FilterType::NoConst;
//! let ignore_const_gaps = false;
//! let ambig_mask = false;
//! ska_array.filter(min_count, filter_ambig_as_missing, &filter, ambig_mask, ignore_const_gaps, update_kmers);
//!
//! // Write out to stdout
//! let mut out_stream = set_ostream(&None);
//! ska_array
//!     .write_fasta(&mut out_stream)
//!     .expect("Couldn't write output fasta");
//! ```
//!
//! To use `ska map` in other rust code, see the example for [`ska_ref`].
//!
//! ## Other APIs
//!
//! It would be possible to make a python API using [`maturin`](https://github.com/PyO3/maturin).
//! Please submit a feature request if you would find this useful.
//!

#![warn(missing_docs)]
#[cfg(any(target_family = "wasm", test))]
use flate2::write::GzEncoder;
#[cfg(any(target_family = "wasm", test))]
use flate2::Compression;
use std::fmt;
#[cfg(any(target_family = "wasm", test))]
use std::io::Write;
#[cfg(not(target_family = "wasm"))]
use std::time::Instant;

use clap::ValueEnum;
extern crate num_cpus;

pub mod merge_ska_dict;
pub mod ska_dict;
#[cfg(not(target_family = "wasm"))]
use crate::merge_ska_dict::build_and_merge;

pub mod ska_ref;
use crate::ska_ref::RefSka;
pub mod merge_ska_array;
#[cfg(not(target_family = "wasm"))]
use crate::merge_ska_array::MergeSkaArray;

#[cfg(not(target_family = "wasm"))]
pub mod generic_modes;
#[cfg(not(target_family = "wasm"))]
use crate::generic_modes::*;

#[cfg(not(target_family = "wasm"))]
pub mod cli;
#[cfg(not(target_family = "wasm"))]
use crate::cli::*;

#[cfg(not(target_family = "wasm"))]
pub mod io_utils;
#[cfg(not(target_family = "wasm"))]
use crate::io_utils::*;

#[cfg(not(target_family = "wasm"))]
pub mod coverage;
#[cfg(not(target_family = "wasm"))]
use crate::coverage::CoverageHistogram;

pub mod cluster;

#[cfg(not(target_family = "wasm"))]
pub mod skalo;
#[cfg(not(target_family = "wasm"))]
use crate::io_utils::load_array;
#[cfg(not(target_family = "wasm"))]
use crate::skalo::utils::Config;

#[cfg(target_family = "wasm")]
use wasm_bindgen::prelude::*;
#[cfg(target_family = "wasm")]
use wasm_bindgen_file_reader::WebSysFile;
#[cfg(target_family = "wasm")]
extern crate console_error_panic_hook;
#[cfg(any(target_family = "wasm", test))]
pub mod wasm;
#[cfg(target_family = "wasm")]
use crate::cluster::{checked_distance_cutoff, cluster_distances_flat};
#[cfg(target_family = "wasm")]
use crate::ska_dict::bit_encoding::UInt;
#[cfg(target_family = "wasm")]
use crate::wasm::ska_align::SkaAlign;
#[cfg(target_family = "wasm")]
use crate::wasm::ska_map::SkaMap;
#[cfg(target_family = "wasm")]
use petgraph::visit::EdgeRef;

/// Possible quality score filters when building with reads
#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
pub enum QualFilter {
    /// Ignore quality scores in reads
    NoFilter,
    /// Filter middle bases below quality threshold
    Middle,
    /// Filter entire k-mer when any base below quality threshold
    Strict,
}

impl fmt::Display for QualFilter {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            Self::NoFilter => write!(f, "No quality filtering"),
            Self::Middle => write!(f, "Middle base quality filtering"),
            Self::Strict => write!(f, "Whole k-mer quality filtering"),
        }
    }
}

/// Quality filtering options for FASTQ files
pub struct QualOpts {
    /// Minimum k-mer count across reads to be added
    pub min_count: u16,
    /// Minimum base quality to be added
    pub min_qual: u8,
    /// [`QualFilter`]: apply quality across whole k-mer, just middle base, or not at all
    pub qual_filter: QualFilter,
}

impl fmt::Display for QualOpts {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "min count: {}; minimum quality {} ({}); filter applied: {}",
            self.min_count,
            self.min_qual,
            (self.min_qual + 33) as char,
            self.qual_filter
        )
    }
}

#[cfg(not(target_family = "wasm"))]
#[doc(hidden)]
pub fn main() {
    let args = cli_args();
    if args.verbose {
        simple_logger::init_with_level(log::Level::Info).unwrap();
    } else {
        simple_logger::init_with_level(log::Level::Warn).unwrap();
    }

    eprintln!("SKA: Split K-mer Analysis (the alignment-free aligner)");
    let start = Instant::now();
    match &args.command {
        Commands::Build {
            seq_files,
            file_list,
            output,
            k,
            single_strand,
            min_count,
            min_qual,
            qual_filter,
            threads,
            proportion_reads,
        } => {
            check_threads(*threads);

            // Read input
            let input_files = get_input_list(file_list, seq_files);
            let rc = !*single_strand;
            // Build, merge
            // check for 64 or 128 bit representation
            let mut quality = QualOpts {
                min_count: 0,
                min_qual: *min_qual,
                qual_filter: *qual_filter,
            };
            if k.le(&31) {
                log::info!("k={}: using 64-bit representation", *k);
                quality.min_count =
                    kmer_min_cutoff::<u64>(min_count, &input_files, k, rc, args.verbose);
                let merged_dict = build_and_merge::<u64>(
                    &input_files,
                    *k,
                    rc,
                    &quality,
                    *threads,
                    *proportion_reads,
                );

                // Save
                save_skf(&merged_dict, output);
            } else {
                log::info!("k={}: using 128-bit representation", *k);
                quality.min_count =
                    kmer_min_cutoff::<u128>(min_count, &input_files, k, rc, args.verbose);
                let merged_dict = build_and_merge::<u128>(
                    &input_files,
                    *k,
                    rc,
                    &quality,
                    *threads,
                    *proportion_reads,
                );

                // Save
                save_skf(&merged_dict, output);
            }
        }
        Commands::Align {
            input,
            output,
            min_freq,
            filter_ambig_as_missing,
            filter,
            ambig_mask,
            no_gap_only_sites,
            threads,
        } => {
            check_threads(*threads);
            if let Ok(mut ska_array) = load_array::<u64>(input, *threads) {
                // In debug mode (cannot be set from CLI, give details)
                log::debug!("{ska_array}");
                align(
                    &mut ska_array,
                    output,
                    filter,
                    *ambig_mask,
                    *no_gap_only_sites,
                    *min_freq,
                    *filter_ambig_as_missing,
                );
            } else if let Ok(mut ska_array) = load_array::<u128>(input, *threads) {
                // In debug mode (cannot be set from CLI, give details)
                log::debug!("{ska_array}");
                align(
                    &mut ska_array,
                    output,
                    filter,
                    *ambig_mask,
                    *no_gap_only_sites,
                    *min_freq,
                    *filter_ambig_as_missing,
                );
            } else {
                panic!("Could not read input file(s): {input:?}");
            }
        }
        Commands::Map {
            reference,
            input,
            output,
            format,
            ambig_mask,
            repeat_mask,
            threads,
        } => {
            check_threads(*threads);

            log::info!("Loading skf as dictionary");
            if let Ok(ska_array) = load_array::<u64>(input, *threads) {
                log::info!(
                    "Making skf of reference k={} rc={}",
                    ska_array.kmer_len(),
                    ska_array.rc()
                );
                let mut ska_ref = RefSka::<u64>::new(
                    ska_array.kmer_len(),
                    reference,
                    ska_array.rc(),
                    *ambig_mask,
                    *repeat_mask,
                );
                map(&ska_array, &mut ska_ref, output, format, *threads);
            } else if let Ok(ska_array) = load_array::<u128>(input, *threads) {
                log::info!(
                    "Making skf of reference k={} rc={}",
                    ska_array.kmer_len(),
                    ska_array.rc()
                );
                let mut ska_ref = RefSka::<u128>::new(
                    ska_array.kmer_len(),
                    reference,
                    ska_array.rc(),
                    *ambig_mask,
                    *repeat_mask,
                );
                map(&ska_array, &mut ska_ref, output, format, *threads);
            } else {
                panic!("Could not read input file(s): {input:?}");
            }
        }
        Commands::Distance {
            skf_file,
            output,
            min_freq,
            allow_ambiguous,
            threads,
            clusters,
            snp_threshold,
        } => {
            check_threads(*threads);
            let filter_ambiguous = !*allow_ambiguous;
            let cluster_threshold = if *clusters {
                Some(*snp_threshold)
            } else {
                None
            };
            if let Ok(mut ska_array) = MergeSkaArray::<u64>::load(skf_file) {
                // In debug mode (cannot be set from CLI, give details)
                log::debug!("{ska_array}");
                distance(
                    &mut ska_array,
                    output,
                    *min_freq,
                    filter_ambiguous,
                    cluster_threshold,
                );
            } else if let Ok(mut ska_array) = MergeSkaArray::<u128>::load(skf_file) {
                // In debug mode (cannot be set from CLI, give details)
                log::debug!("{ska_array}");
                distance(
                    &mut ska_array,
                    output,
                    *min_freq,
                    filter_ambiguous,
                    cluster_threshold,
                );
            } else {
                panic!("Could not read input file(s): {skf_file}");
            }
        }
        Commands::Merge { skf_files, output } => {
            if skf_files.len() < 2 {
                panic!("Need at least two files to merge");
            }

            log::info!("Loading first alignment");
            if let Ok(first_array) = MergeSkaArray::<u64>::load(&skf_files[0]) {
                merge(&first_array, &skf_files[1..], output);
            } else if let Ok(first_array) = MergeSkaArray::<u128>::load(&skf_files[0]) {
                merge(&first_array, &skf_files[1..], output);
            } else {
                panic!("Could not read input file: {}", skf_files[0]);
            }
        }
        Commands::Delete {
            skf_file,
            output,
            file_list,
            names,
        } => {
            let input_files = get_input_list(file_list, names);
            let input_names: Vec<&str> = input_files.iter().map(|t| &*t.0).collect();
            let output_file = output.clone().unwrap_or(skf_file.to_string());
            log::info!("Loading skf file");
            if let Ok(mut ska_array) = MergeSkaArray::<u64>::load(skf_file) {
                delete(&mut ska_array, &input_names, &output_file);
            } else if let Ok(mut ska_array) = MergeSkaArray::<u128>::load(skf_file) {
                delete(&mut ska_array, &input_names, &output_file);
            } else {
                panic!("Could not read input file: {skf_file}");
            }
        }
        Commands::Weed {
            skf_file,
            weed_file,
            output,
            reverse,
            min_freq,
            filter_ambig_as_missing,
            ambig_mask,
            no_gap_only_sites,
            filter,
        } => {
            log::info!("Loading skf file");
            if let Ok(mut ska_array) = MergeSkaArray::<u64>::load(skf_file) {
                weed(
                    &mut ska_array,
                    weed_file,
                    *reverse,
                    *min_freq,
                    *filter_ambig_as_missing,
                    filter,
                    *ambig_mask,
                    *no_gap_only_sites,
                    if output.is_none() {
                        skf_file
                    } else {
                        output.as_ref().unwrap().as_str()
                    },
                );
            } else if let Ok(mut ska_array) = MergeSkaArray::<u128>::load(skf_file) {
                weed(
                    &mut ska_array,
                    weed_file,
                    *reverse,
                    *min_freq,
                    *filter_ambig_as_missing,
                    filter,
                    *ambig_mask,
                    *no_gap_only_sites,
                    if output.is_none() {
                        skf_file
                    } else {
                        output.as_ref().unwrap().as_str()
                    },
                );
            } else {
                panic!("Could not read input file: {skf_file}");
            }
        }
        Commands::Nk {
            skf_file,
            full_info,
        } => {
            if let Ok(ska_array) = MergeSkaArray::<u64>::load(skf_file) {
                println!("{ska_array}");
                if *full_info {
                    log::info!("Printing full info");
                    println!("{ska_array:?}");
                }
            } else if let Ok(ska_array) = MergeSkaArray::<u128>::load(skf_file) {
                println!("{ska_array}");
                if *full_info {
                    log::info!("Printing full info");
                    println!("{ska_array:?}");
                }
            } else {
                panic!("Could not read input file: {skf_file}");
            }
        }
        Commands::Cov {
            fastq_fwd,
            fastq_rev,
            k,
            single_strand,
        } => {
            // Build, merge
            let rc = !*single_strand;
            let cutoff;
            if *k <= 31 {
                log::info!("k={}: using 64-bit representation", *k);
                let mut cov =
                    CoverageHistogram::<u64>::new(fastq_fwd, fastq_rev, *k, rc, args.verbose);
                cutoff = cov.fit_histogram().expect("Couldn't fit coverage model");
                cov.plot_hist();
            } else {
                log::info!("k={}: using 128-bit representation", *k);
                let mut cov =
                    CoverageHistogram::<u128>::new(fastq_fwd, fastq_rev, *k, rc, args.verbose);
                cutoff = cov.fit_histogram().expect("Couldn't fit coverage model");
                cov.plot_hist();
            }
            eprintln!("Estimated cutoff\t{cutoff}");
        }
        Commands::Lo {
            input_skf,
            reference,
            output,
            missing,
            depth,
            indel_kmers,
            threads,
        } => {
            let config = Config {
                input_file: input_skf.clone(),
                output_name: output.clone(),
                max_missing: *missing,
                max_depth: *depth,
                max_indel_kmers: *indel_kmers,
                nb_threads: *threads,
                reference_genome: reference.clone(),
            };

            if let Ok(ska_array) = load_array::<u64>(std::slice::from_ref(input_skf), *threads) {
                log::info!("Reading file {input_skf}");
                log::info!("Using 64-bit representation");
                skalo(ska_array, config);
            } else if let Ok(ska_array) =
                load_array::<u128>(std::slice::from_ref(input_skf), *threads)
            {
                log::info!("Reading file {input_skf}");
                log::info!("Using 128-bit representation");
                skalo(ska_array, config);
            } else {
                panic!("Could not read input file(s): {input_skf:?}");
            }
        }
    }
    let end = Instant::now();

    eprintln!("SKA done in {}s", end.duration_since(start).as_secs());
    eprintln!("⬛⬜⬛⬜⬛⬜⬛");
    eprintln!("⬜⬛⬜⬛⬜⬛⬜");
    log::info!("Complete");
}

// WASM implementation
#[cfg(target_family = "wasm")]
#[doc(hidden)]
pub fn main() {
    panic!("You've compiled Ska2 for WebAssembly support, you cannot use it as a normal binary anymore!");
}

#[cfg(target_family = "wasm")]
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(js_namespace = console)]
    fn log(s: &str);

    #[wasm_bindgen(js_name = postMessage)]
    fn post_message(data: &JsValue);
}

/// This allows to send a JS object with data included as message, altogether, along the post_message method previously defined.
#[cfg(target_family = "wasm")]
pub(crate) fn post_alignment_progress(
    run_id: u32,
    stage: &str,
    sample_index: Option<usize>,
    sample_total: Option<usize>,
    sample_name: Option<&str>,
) {
    let progress = js_sys::Object::new();
    let _ = js_sys::Reflect::set(
        &progress,
        &"runId".into(),
        &JsValue::from_f64(run_id as f64),
    );
    let _ = js_sys::Reflect::set(&progress, &"stage".into(), &JsValue::from_str(stage));
    let _ = js_sys::Reflect::set(
        &progress,
        &"sampleIndex".into(),
        &sample_index.map_or(JsValue::NULL, |index| JsValue::from_f64(index as f64)),
    );
    let _ = js_sys::Reflect::set(
        &progress,
        &"sampleTotal".into(),
        &sample_total.map_or(JsValue::NULL, |total| JsValue::from_f64(total as f64)),
    );
    let _ = js_sys::Reflect::set(
        &progress,
        &"sampleName".into(),
        &sample_name.map_or(JsValue::NULL, JsValue::from_str),
    );

    let message = js_sys::Object::new();
    let _ = js_sys::Reflect::set(&message, &"alignmentProgress".into(), &progress);
    post_message(&message.into());
}

#[cfg(target_family = "wasm")]
#[wasm_bindgen]
/// Function that allows to propagate panic error messages when compiling to wasm, see https://github.com/rustwasm/console_error_panic_hook
pub fn init_panic_hook() {
    console_error_panic_hook::set_once();
}

/// Logging wrapper function for the WebAssembly version
#[cfg(target_family = "wasm")]
pub fn logw(text: &str, typ: Option<&str>) {
    if let Some(thetyp) = typ {
        log((String::from("ska.rust::") + thetyp + "::" + text).as_str());
    } else {
        log(text);
    }
}

#[cfg(target_family = "wasm")]
#[wasm_bindgen]
/// Struct to interact with JS when working with WebAssembly
pub struct SkaData {
    k: usize,
    rc: bool,
    reference64: Option<RefSka<u64>>,
    reference128: Option<RefSka<u128>>,
    reference_string: Vec<String>,
    mapped64: Option<Vec<SkaMap<u64>>>,
    mapped128: Option<Vec<SkaMap<u128>>>,
    sample_names: Vec<String>,
}

#[cfg(target_family = "wasm")]
#[wasm_bindgen]
impl SkaData {
    /// Constructor of the SkaData struct
    pub fn new(
        ref_file: web_sys::File,
        k: u32,
        rc: bool,
        ambig_mask: bool,
        repeat_mask: bool,
    ) -> Self {
        let k = k as usize;

        if cfg!(debug_assertions) {
            init_panic_hook();
        }

        let mut wf = WebSysFile::new(ref_file);
        if k < 32 {
            let reference = RefSka::<u64>::new(k, &mut wf, rc, ambig_mask, repeat_mask);
            logw(
                &format!("Indexed reference: {} split k-mers", reference.ksize()),
                None,
            );
            let reference_string = reconstruct_sequence(&reference);

            Self {
                k,
                rc,
                reference64: Some(reference),
                reference128: None,
                reference_string,
                mapped64: Some(Vec::new()),
                mapped128: None,
                sample_names: Vec::new(),
            }
        } else if k < 64 {
            let reference = RefSka::<u128>::new(k, &mut wf, rc, ambig_mask, repeat_mask);
            logw(
                &format!("Indexed reference: {} split k-mers", reference.ksize()),
                None,
            );
            let reference_string = reconstruct_sequence(&reference);

            Self {
                k,
                rc,
                reference64: None,
                reference128: Some(reference),
                reference_string,
                mapped64: None,
                mapped128: Some(Vec::new()),
                sample_names: Vec::new(),
            }
        } else {
            panic!("k values larger than 64 not supported");
        }
    }

    /// Mapping function.
    pub fn map(
        &mut self,
        input_file: web_sys::File,
        rev_reads: Option<web_sys::File>,
        proportion_reads: Option<f64>,
        min_count: u16,
        min_qual: u8,
        qual_filter: u8,
        sample_name: String,
    ) -> String {
        if rev_reads.is_some() {
            logw("Detected paired fastq input files", None);
        }
        logw("Mapping reads to reference", None);
        let input_file_ref = input_file;
        if self.k < 32 {
            if let Some(second_file) = rev_reads {
                self.mapped64.as_mut().unwrap().push(SkaMap::<u64>::new(
                    self.reference64.as_ref().unwrap(),
                    &input_file_ref,
                    Some(&second_file),
                    proportion_reads,
                    self.rc,
                    min_count,
                    min_qual,
                    if qual_filter == 0 {
                        QualFilter::NoFilter
                    } else if qual_filter == 1 {
                        QualFilter::Middle
                    } else {
                        QualFilter::Strict
                    },
                ));
            } else {
                self.mapped64.as_mut().unwrap().push(SkaMap::<u64>::new(
                    self.reference64.as_ref().unwrap(),
                    &input_file_ref,
                    None,
                    proportion_reads,
                    self.rc,
                    min_count,
                    min_qual,
                    if qual_filter == 0 {
                        QualFilter::NoFilter
                    } else if qual_filter == 1 {
                        QualFilter::Middle
                    } else {
                        QualFilter::Strict
                    },
                ));
            };
        } else if let Some(second_file) = rev_reads {
            self.mapped128.as_mut().unwrap().push(SkaMap::<u128>::new(
                self.reference128.as_ref().unwrap(),
                &input_file_ref,
                Some(&second_file),
                proportion_reads,
                self.rc,
                min_count,
                min_qual,
                if qual_filter == 0 {
                    QualFilter::NoFilter
                } else if qual_filter == 1 {
                    QualFilter::Middle
                } else {
                    QualFilter::Strict
                },
            ));
        } else {
            self.mapped128.as_mut().unwrap().push(SkaMap::<u128>::new(
                self.reference128.as_ref().unwrap(),
                &input_file_ref,
                None,
                proportion_reads,
                self.rc,
                min_count,
                min_qual,
                if qual_filter == 0 {
                    QualFilter::NoFilter
                } else if qual_filter == 1 {
                    QualFilter::Middle
                } else {
                    QualFilter::Strict
                },
            ));
        }

        self.sample_names.push(sample_name);

        logw("Reads mapped successfully!", None);

        let mut results = json::JsonValue::new_array();

        results["Mapped sequences"] = json::JsonValue::new_array();
        let whole_mapping = if self.k < 32 {
            self.reference64.as_ref().unwrap().pseudoalignment(
                self.mapped64.as_ref().unwrap()[self.mapped64.as_ref().unwrap().len() - 1]
                    .get_mapped_bases(),
            )[0]
            .clone()
        } else {
            self.reference128.as_ref().unwrap().pseudoalignment(
                self.mapped128.as_ref().unwrap()[self.mapped128.as_ref().unwrap().len() - 1]
                    .get_mapped_bases(),
            )[0]
            .clone()
        };

        let mut current_length = 0;
        for chr in 0..self.reference_string.len() {
            let chr_mapping: String = whole_mapping
                [current_length..current_length + self.reference_string[chr].len()]
                .to_string();
            let _ = results["Mapped sequences"].push(chr_mapping);
            current_length += self.reference_string[chr].len();
        }

        if self.k < 32 {
            results["Number of variants"] = self.mapped64.as_ref().unwrap()
                [self.mapped64.as_ref().unwrap().len() - 1]
                .get_mapped_bases()
                .len()
                .into();
        } else {
            results["Number of variants"] = self.mapped128.as_ref().unwrap()
                [self.mapped128.as_ref().unwrap().len() - 1]
                .get_mapped_bases()
                .len()
                .into();
        }
        // The result is a concatenated string of all the mapped bases

        let mut count_mapped_bases = 0;
        let mut count_total_bases = 0;
        for base in whole_mapping.chars() {
            if base != '-' {
                count_mapped_bases += 1;
                count_total_bases += 1;
            } else {
                count_total_bases += 1;
            }
        }

        results["Coverage"] = (count_mapped_bases as f64 / count_total_bases as f64).into();
        results["VCF"] = self.current_vcf_text().into();

        logw("Results computed successfully!", None);

        results.dump()
    }

    fn current_vcf_text(&self) -> String {
        if self.k < 32 {
            let sample_maps: Vec<Vec<_>> = self
                .mapped64
                .as_ref()
                .unwrap()
                .iter()
                .map(|mapped| mapped.get_mapped_bases().clone())
                .collect();
            self.reference64
                .as_ref()
                .unwrap()
                .write_vcf_text_from_sparse_maps(&self.sample_names, &sample_maps)
        } else {
            let sample_maps: Vec<Vec<_>> = self
                .mapped128
                .as_ref()
                .unwrap()
                .iter()
                .map(|mapped| mapped.get_mapped_bases().clone())
                .collect();
            self.reference128
                .as_ref()
                .unwrap()
                .write_vcf_text_from_sparse_maps(&self.sample_names, &sample_maps)
        }
    }

    /// Retrieves the reference sequence
    pub fn get_reference(&self) -> String {
        self.reference_string.join("\n")
    }
}

#[cfg(target_family = "wasm")]
/// Reconstructs a sequence beginning from a reference
pub fn reconstruct_sequence<IntT: for<'a> UInt<'a>>(reference: &RefSka<IntT>) -> Vec<String> {
    let sequence_u8 = reference.get_seq();
    let mut sequence_string = Vec::new();

    for sequence in sequence_u8 {
        let mut current_seq = "".to_string();
        for base in sequence {
            if *base != 10 {
                current_seq.push(*base as char);
            }
        }
        sequence_string.push(current_seq);
    }

    sequence_string
}

#[cfg(target_family = "wasm")]
#[wasm_bindgen]
/// AlignData struct for doing alignment while in WebAssembly
pub struct AlignData {
    run_id: u32,
    k: usize,
    alignment64: Option<SkaAlign<u64>>,
    alignment128: Option<SkaAlign<u128>>,
    file_names: Vec<String>,
    sample_total: usize,
}

#[cfg(target_family = "wasm")]
impl AlignData {
    fn tree_result(&self, on_stage: &dyn Fn(&str)) -> Result<json::JsonValue, JsValue> {
        if self.file_names.len() != self.sample_total {
            return Err(JsValue::from_str(
                "Not all planned samples have been added.",
            ));
        }
        let size = if self.k < 32 {
            self.alignment64
                .as_ref()
                .ok_or_else(|| JsValue::from_str("No 64-bit alignment is available."))?
                .get_size()
        } else {
            self.alignment128
                .as_ref()
                .ok_or_else(|| JsValue::from_str("No 128-bit alignment is available."))?
                .get_size()
        };
        if size != self.file_names.len() {
            return Err(JsValue::from_str(
                "Distance sample count does not match the sample names.",
            ));
        }
        let tree_tip_names = if size >= 3 {
            crate::wasm::ska_align::make_tree_tip_names(&self.file_names)
        } else {
            Vec::new()
        };
        let newick = if size < 3 {
            String::new()
        } else {
            let result = if self.k < 32 {
                self.alignment64
                    .as_ref()
                    .ok_or_else(|| JsValue::from_str("No 64-bit alignment is available."))?
                    .build_tree(&tree_tip_names, on_stage)
            } else {
                self.alignment128
                    .as_ref()
                    .ok_or_else(|| JsValue::from_str("No 128-bit alignment is available."))?
                    .build_tree(&tree_tip_names, on_stage)
            };
            result.map_err(|error| JsValue::from_str(&error))?
        };
        let mut results = json::JsonValue::new_object();
        results["newick"] = newick.into();
        results["names"] = json::JsonValue::new_array();
        results["treeTipNames"] = json::JsonValue::new_array();
        for name in &self.file_names {
            results["names"]
                .push(name.clone())
                .map_err(|error| JsValue::from_str(&error.to_string()))?;
        }
        for name in tree_tip_names {
            results["treeTipNames"]
                .push(name)
                .map_err(|error| JsValue::from_str(&error.to_string()))?;
        }
        Ok(results)
    }
}

/// Get the distances and compress them to gz directly.
#[cfg(any(target_family = "wasm", test))]
fn distances_csv_gzip(names: &[String], flat_distances: &[u32]) -> std::io::Result<Vec<u8>> {
    if names.len() < 3 {
        return Ok(Vec::new());
    }

    let n = names.len();
    let expected = crate::wasm::alignment_matrix::strict_triangle_len(n).ok_or_else(|| {
        std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "Sample count exceeds the distance index range.",
        )
    })?;
    if flat_distances.len() != expected {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "Distance count does not match the sample count.",
        ));
    }

    let mut encoder = GzEncoder::new(Vec::new(), Compression::fast());
    encoder.write_all(b",")?;
    for (index, name) in names.iter().enumerate() {
        if index > 0 {
            encoder.write_all(b",")?;
        }
        encoder.write_all(name.as_bytes())?;
    }
    encoder.write_all(b"\n")?;

    for (i, name) in names.iter().enumerate() {
        encoder.write_all(name.as_bytes())?;
        for j in 0..n {
            encoder.write_all(b",")?;
            if i == j {
                encoder.write_all(b"0")?;
            } else {
                let (row, column) = if i < j { (j, i) } else { (i, j) };
                let flat_idx = crate::wasm::alignment_matrix::strict_triangle_index(row, column)
                    .ok_or_else(|| {
                        std::io::Error::new(
                            std::io::ErrorKind::InvalidData,
                            "Invalid strict-triangle distance index.",
                        )
                    })?;
                let distance = flat_distances[flat_idx];
                write!(encoder, "{distance}")?;
            }
        }
        encoder.write_all(b"\n")?;
    }

    encoder.finish()
}

#[cfg(target_family = "wasm")]
#[wasm_bindgen]
/// Extracts one sample's split-kmer dictionary into compact, transferable words.
#[allow(clippy::too_many_arguments)]
pub fn extract_alignment_sample(
    file1: web_sys::File,
    file2: Option<web_sys::File>,
    k: u32,
    sample_index: u32,
    sample_name: String,
    rc: bool,
    proportion_reads: Option<f64>,
    min_count: u16,
    min_qual: u8,
    qual_filter: u8,
) -> Result<js_sys::Uint32Array, JsValue> {
    init_panic_hook();
    let qual_filter = match qual_filter {
        0 => QualFilter::NoFilter,
        1 => QualFilter::Middle,
        _ => QualFilter::Strict,
    };
    let qual = QualOpts {
        min_count,
        min_qual,
        qual_filter,
    };

    logw(
        &format!("Extracting split k-mers for sample {}", sample_name),
        Some("info"),
    );
    let words = if k < 32 {
        let dictionary = crate::ska_dict::SkaDict::<u64>::new(
            k as usize,
            sample_index as usize,
            (&file1, file2.as_ref()),
            &sample_name,
            rc,
            &qual,
            proportion_reads,
        );
        let mut words = Vec::with_capacity(dictionary.ksize() * 2);
        for (kmer, base) in dictionary.into_kmers() {
            words.extend(crate::wasm::alignment_matrix::pack_entry_u64(kmer, base));
        }
        words
    } else {
        let dictionary = crate::ska_dict::SkaDict::<u128>::new(
            k as usize,
            sample_index as usize,
            (&file1, file2.as_ref()),
            &sample_name,
            rc,
            &qual,
            proportion_reads,
        );
        let mut words = Vec::with_capacity(dictionary.ksize() * 4);
        for (kmer, base) in dictionary.into_kmers() {
            words.extend(crate::wasm::alignment_matrix::pack_entry_u128(kmer, base));
        }
        words
    };

    logw(
        &format!(
            "Extracted {} split k-mers for sample {}",
            words.len() / if k < 32 { 2 } else { 4 },
            sample_name
        ),
        Some("info"),
    );
    Ok(js_sys::Uint32Array::from(words.as_slice()))
}

#[cfg(target_family = "wasm")]
struct EncodedAlignmentSequence {
    name: String,
    lo: Vec<u64>,
    hi: Vec<u64>,
    valid: Vec<u64>,
}

#[cfg(target_family = "wasm")]
fn concrete_base_bits(base: u8) -> Option<(bool, bool)> {
    crate::wasm::alignment_matrix::concrete_base_code(base)
        .map(|code| (code & 1 != 0, code & 2 != 0))
}

#[cfg(target_family = "wasm")]
fn encode_alignment_sequence(
    name: String,
    sequence: &[u8],
    n_words: usize,
) -> Result<EncodedAlignmentSequence, JsValue> {
    let mut encoded = EncodedAlignmentSequence {
        name: name.clone(),
        lo: vec![0; n_words],
        hi: vec![0; n_words],
        valid: vec![0; n_words],
    };

    for (pos, base) in sequence.iter().enumerate() {
        if !base.is_ascii() {
            return Err(JsValue::from_str(&format!(
                "Alignment sequence for {name} contains non-ASCII characters."
            )));
        }

        if let Some((lo_bit, hi_bit)) = concrete_base_bits(*base) {
            let word = pos / 64;
            let mask = 1u64 << (pos % 64);
            encoded.valid[word] |= mask;
            if lo_bit {
                encoded.lo[word] |= mask;
            }
            if hi_bit {
                encoded.hi[word] |= mask;
            }
        }
    }

    Ok(encoded)
}

#[cfg(any(target_family = "wasm", test))]
fn add_imported_mismatches(distance: u32, increment: u32) -> (u32, bool) {
    match distance.checked_add(increment) {
        Some(total) => (total, false),
        None => (u32::MAX, true),
    }
}

#[cfg(target_family = "wasm")]
fn packed_snp_distance(a: &EncodedAlignmentSequence, b: &EncodedAlignmentSequence) -> (u32, bool) {
    let mut distance = 0u32;
    let mut capped = false;

    for word_idx in 0..a.valid.len() {
        let comparable = a.valid[word_idx] & b.valid[word_idx];
        let lo_diff = a.lo[word_idx] ^ b.lo[word_idx];
        let hi_diff = a.hi[word_idx] ^ b.hi[word_idx];
        let snp_bits = (lo_diff | hi_diff) & comparable;
        let (new_distance, overflowed) = add_imported_mismatches(distance, snp_bits.count_ones());
        distance = new_distance;
        capped |= overflowed;
    }

    (distance, capped)
}

#[cfg(target_family = "wasm")]
fn parse_fasta_alignment(text: &str) -> Result<Vec<(String, Vec<u8>)>, JsValue> {
    let mut records: Vec<(String, Vec<u8>)> = Vec::new();
    let mut current_name: Option<String> = None;
    let mut current_sequence: Vec<u8> = Vec::new();
    let mut seen_names = std::collections::HashSet::new();
    for raw_line in text.lines() {
        let line = raw_line.trim();
        if line.is_empty() {
            continue;
        }
        if let Some(header) = line.strip_prefix('>') {
            if let Some(name) = current_name.take() {
                if current_sequence.is_empty() {
                    return Err(JsValue::from_str(&format!(
                        "Alignment record {name} has no sequence."
                    )));
                }
                records.push((name, std::mem::take(&mut current_sequence)));
            }

            let name = header
                .split_whitespace()
                .next()
                .ok_or_else(|| JsValue::from_str("FASTA headers must include a sample ID."))?
                .to_string();
            if name.is_empty() {
                return Err(JsValue::from_str("FASTA headers must include a sample ID."));
            }
            if !seen_names.insert(name.clone()) {
                return Err(JsValue::from_str(&format!(
                    "Duplicate sample ID in alignment: {name}."
                )));
            }
            current_name = Some(name);
        } else {
            if current_name.is_none() {
                return Err(JsValue::from_str(
                    "Alignment files must be FASTA-formatted and start with a '>' header.",
                ));
            }
            for base in line.bytes() {
                if !base.is_ascii() {
                    let name = current_name.as_deref().unwrap_or("unknown sample");
                    return Err(JsValue::from_str(&format!(
                        "Alignment sequence for {name} contains non-ASCII characters."
                    )));
                }
                if !base.is_ascii_whitespace() {
                    current_sequence.push(base.to_ascii_uppercase());
                }
            }
        }
    }

    if let Some(name) = current_name.take() {
        if current_sequence.is_empty() {
            return Err(JsValue::from_str(&format!(
                "Alignment record {name} has no sequence."
            )));
        }
        records.push((name, current_sequence));
    }

    if records.len() < 2 {
        return Err(JsValue::from_str(
            "Transmission clustering requires at least two aligned sequences.",
        ));
    }

    let alignment_len = records[0].1.len();
    if alignment_len == 0 {
        return Err(JsValue::from_str("Alignment sequences must not be empty."));
    }

    for (name, sequence) in &records {
        if sequence.len() != alignment_len {
            return Err(JsValue::from_str(&format!(
                "Alignment sequence {name} has length {}, expected {alignment_len}.",
                sequence.len()
            )));
        }
    }

    Ok(records)
}

/// Return a normalised pair key, the mate number, and a sample name for a FASTQ
/// filename with a recognised mate suffix (`_1`, `_2`, `_R1`, or `_R2`).
#[cfg(any(target_family = "wasm", test))]
fn fastq_pair_info(filename: &str) -> Option<(String, u8, String)> {
    let without_gzip = filename.strip_suffix(".gz").unwrap_or(filename);
    let extension_start = without_gzip.rfind('.')?;
    let extension = &without_gzip[extension_start + 1..];
    if !["fastq", "fq", "fnq"].contains(&extension) {
        return None;
    }

    let stem = &without_gzip[..extension_start];
    let bytes = stem.as_bytes();
    let is_delimiter = |byte| matches!(byte, b'_' | b'.' | b'-');

    // Prefer the explicit R1/R2 marker, which may be followed by an Illumina
    // chunk suffix such as `_001`.
    let mut marker = None;
    for index in 1..bytes.len().saturating_sub(1) {
        if bytes[index] == b'R'
            && is_delimiter(bytes[index - 1])
            && matches!(bytes[index + 1], b'1' | b'2')
        {
            let suffix = &stem[index + 2..];
            let valid_suffix = suffix.is_empty()
                || suffix.strip_prefix('_').is_some_and(|digits| {
                    !digits.is_empty() && digits.bytes().all(|b| b.is_ascii_digit())
                });
            if valid_suffix {
                marker = Some((index - 1, index + 1, index + 2, bytes[index + 1] - b'0'));
            }
        }
    }

    // Also accept the common terminal `_1` / `_2` and `.1` / `.2` forms.
    if marker.is_none() && bytes.len() >= 2 {
        let digit_index = bytes.len() - 1;
        if matches!(bytes[digit_index], b'1' | b'2') && is_delimiter(bytes[digit_index - 1]) {
            marker = Some((
                digit_index - 1,
                digit_index,
                digit_index + 1,
                bytes[digit_index] - b'0',
            ));
        }
    }

    let (marker_start, mate_index, marker_end, mate) = marker?;
    let mut key = bytes.to_vec();
    key[mate_index] = b'X';
    let pair_key = String::from_utf8(key).ok()?;
    let sample_stem = format!("{}{}", &stem[..marker_start], &stem[marker_end..]);
    let sample_name = format!("{sample_stem}{}", &filename[stem.len()..]);

    Some((pair_key, mate, sample_name))
}

#[cfg(any(target_family = "wasm", test))]
fn fastq_files_are_mates(first: &str, second: &str) -> bool {
    matches!(
        (fastq_pair_info(first), fastq_pair_info(second)),
        (Some((key1, mate1, _)), Some((key2, mate2, _)))
            if key1 == key2 && mate1 != mate2
    )
}

#[cfg(any(target_family = "wasm", test))]
#[derive(Debug, PartialEq, Eq)]
struct AlignmentSample {
    name: String,
    first_index: usize,
    second_index: Option<usize>,
}

#[cfg(any(target_family = "wasm", test))]
fn is_fastq_filename(filename: &str) -> bool {
    let mut file_type = filename.rsplit('.').next().unwrap_or("");
    if file_type == "gz" {
        file_type = filename.rsplit('.').nth(1).unwrap_or("");
    }
    ["fnq", "fq", "fastq"].contains(&file_type)
}

#[cfg(any(target_family = "wasm", test))]
fn alignment_sample_for_read(names: &[String], index: usize) -> AlignmentSample {
    AlignmentSample {
        name: names[index].clone(),
        first_index: index,
        second_index: None,
    }
}

/// Group the selected inputs in the same order as the alignment workflow:
/// non-FASTQ inputs first, then paired or individual FASTQ samples.
#[cfg(any(target_family = "wasm", test))]
fn prepare_alignment_samples(names: &[String]) -> Vec<AlignmentSample> {
    let mut samples = Vec::with_capacity(names.len());
    let mut fastq_files = Vec::new();

    for (index, filename) in names.iter().enumerate() {
        if is_fastq_filename(filename) {
            fastq_files.push(index);
        } else {
            samples.push(AlignmentSample {
                name: filename.clone(),
                first_index: index,
                second_index: None,
            });
        }
    }

    match fastq_files.len() {
        0 => {}
        1 => {
            #[cfg(target_family = "wasm")]
            logw("Among the uploaded files, only one fastq file has been uploaded. We'll assume it contains all reads needed for this.", None);
            samples.push(alignment_sample_for_read(names, fastq_files[0]));
        }
        2 => {
            let first = fastq_files[0];
            let second = fastq_files[1];
            if fastq_files_are_mates(&names[first], &names[second]) {
                let sample_name = fastq_pair_info(&names[first])
                    .expect("paired FASTQ filename should have mate information")
                    .2;
                samples.push(AlignmentSample {
                    name: sample_name,
                    first_index: first,
                    second_index: Some(second),
                });
            } else {
                #[cfg(target_family = "wasm")]
                logw(
                    "Two read files uploaded, but found to not come from the same sample.",
                    None,
                );
                samples.push(alignment_sample_for_read(names, first));
                samples.push(alignment_sample_for_read(names, second));
            }
        }
        _ => {
            while let Some(first) = fastq_files.pop() {
                let mate_position = fastq_files
                    .iter()
                    .position(|&second| fastq_files_are_mates(&names[first], &names[second]));

                if let Some(position) = mate_position {
                    let second = fastq_files.remove(position);
                    #[cfg(target_family = "wasm")]
                    logw(
                        &format!("Matched FASTQ pair: {} and {}", names[first], names[second]),
                        Some("info"),
                    );
                    let sample_name = fastq_pair_info(&names[first])
                        .expect("paired FASTQ filename should have mate information")
                        .2;
                    samples.push(AlignmentSample {
                        name: sample_name,
                        first_index: first,
                        second_index: Some(second),
                    });
                } else {
                    #[cfg(target_family = "wasm")]
                    logw(
                        &format!(
                            "No pair found among input files for {}. Adding individually.",
                            names[first]
                        ),
                        None,
                    );
                    samples.push(alignment_sample_for_read(names, first));
                }
            }
        }
    }

    samples
}

#[cfg(target_family = "wasm")]
#[wasm_bindgen]
impl AlignData {
    /// Constructor of the AlignData struct
    pub fn new(k: u32, run_id: u32) -> Self {
        init_panic_hook();
        let k = k as usize;

        if k < 32 {
            Self {
                run_id,
                k,
                alignment64: Some(SkaAlign::<u64>::new()),
                alignment128: None,
                file_names: Vec::new(),
                sample_total: 0,
            }
        } else if k < 64 {
            Self {
                run_id,
                k,
                alignment64: None,
                alignment128: Some(SkaAlign::<u128>::new()),
                file_names: Vec::new(),
                sample_total: 0,
            }
        } else {
            panic!("k values larger than 64 are not supported.");
        }
    }

    /// Builds clustering-ready alignment data from an already aligned FASTA file.
    pub fn from_alignment_text(text: String) -> Result<AlignData, JsValue> {
        let records = parse_fasta_alignment(&text)?;
        let alignment_len = records[0].1.len();
        let n_words = alignment_len
            .checked_add(63)
            .ok_or_else(|| JsValue::from_str("Alignment length is too large."))?
            / 64;
        let mut encoded = Vec::with_capacity(records.len());

        for (name, sequence) in records {
            encoded.push(encode_alignment_sequence(name, &sequence, n_words)?);
        }

        let mut file_names = Vec::with_capacity(encoded.len());
        for sequence in &encoded {
            file_names.push(sequence.name.clone());
        }

        let n = encoded.len();
        let pair_count = crate::wasm::alignment_matrix::strict_triangle_len(n)
            .ok_or_else(|| JsValue::from_str("Imported alignment has too many samples."))?;
        let mut distances = Vec::new();
        distances.try_reserve_exact(pair_count).map_err(|error| {
            JsValue::from_str(&format!("Could not allocate distances: {error}"))
        })?;
        let mut distances_capped = false;
        for row in 1..n {
            for column in 0..row {
                let (distance, capped) = packed_snp_distance(&encoded[row], &encoded[column]);
                distances.push(distance);
                distances_capped |= capped;
            }
        }
        let imported_alignment = SkaAlign::<u64>::from_imported_distances(
            file_names.clone(),
            distances,
            distances_capped,
        )
        .map_err(|error| JsValue::from_str(&error))?;

        Ok(Self {
            run_id: 0,
            k: 31,
            alignment64: Some(imported_alignment),
            alignment128: None,
            file_names,
            sample_total: n,
        })
    }

    /// Plan an initial batch or append, validating before modifying the dataset.
    pub fn plan_samples(
        &mut self,
        input_names: Vec<String>,
        run_id: u32,
    ) -> Result<String, JsValue> {
        post_alignment_progress(run_id, "preparing", None, None, None);
        let frozen = self
            .alignment64
            .as_ref()
            .map(|a| a.is_sealed())
            .or_else(|| self.alignment128.as_ref().map(|a| a.is_sealed()))
            .ok_or_else(|| JsValue::from_str("No extracted alignment dataset is available."))?;
        if frozen {
            return Err(JsValue::from_str("The alignment is closed to new samples."));
        }
        let samples = prepare_alignment_samples(&input_names);
        if samples.is_empty() {
            return Err(JsValue::from_str("No supported input samples were found."));
        }
        let total = self
            .file_names
            .len()
            .checked_add(samples.len())
            .filter(|total| *total <= crate::wasm::alignment_matrix::MAX_SAMPLES)
            .ok_or_else(|| {
                JsValue::from_str(&format!(
                    "Sparse alignment supports at most {} samples.",
                    crate::wasm::alignment_matrix::MAX_SAMPLES
                ))
            })?;
        let mut names: std::collections::HashSet<&str> =
            self.file_names.iter().map(String::as_str).collect();
        for sample in &samples {
            if !names.insert(&sample.name) {
                return Err(JsValue::from_str(&format!("Duplicate sample name: {}. Upload each sample only once; paired mates must arrive together.", sample.name)));
            }
        }
        let mut plan = json::JsonValue::new_array();
        for (offset, sample) in samples.iter().enumerate() {
            let mut item = json::JsonValue::new_object();
            item["sampleIndex"] = (self.file_names.len() + offset).into();
            item["name"] = sample.name.clone().into();
            item["firstFileIndex"] = sample.first_index.into();
            item["secondFileIndex"] = sample
                .second_index
                .map_or(json::JsonValue::Null, |index| index.into());
            let _ = plan.push(item);
        }
        if let Some(alignment) = self.alignment64.as_mut() {
            alignment
                .set_expected_samples(total)
                .map_err(|error| JsValue::from_str(&error))?;
        }
        if let Some(alignment) = self.alignment128.as_mut() {
            alignment
                .set_expected_samples(total)
                .map_err(|error| JsValue::from_str(&error))?;
        }
        self.run_id = run_id;
        self.sample_total = total;
        logw(
            &format!(
                "Prepared {} additional sample(s); {} samples in total",
                samples.len(),
                total
            ),
            Some("info"),
        );
        Ok(plan.dump())
    }

    /// Adds one sample's compact transferred split-kmer entries in sample order.
    pub fn add_sample_words(
        &mut self,
        sample_index: u32,
        sample_name: String,
        words: Vec<u32>,
    ) -> Result<(), JsValue> {
        let expected_index = self.file_names.len();
        if expected_index >= self.sample_total {
            return Err(JsValue::from_str(
                "The sample is outside the validated batch plan.",
            ));
        }
        if sample_index as usize != expected_index {
            return Err(JsValue::from_str(&format!(
                "Alignment sample order mismatch: received {}, expected {}.",
                sample_index, expected_index
            )));
        }

        let words_per_entry = if self.k < 32 { 2 } else { 4 };
        if words.len() % words_per_entry != 0 {
            return Err(JsValue::from_str(&format!(
                "Packed sample {} has {} words; expected a multiple of {}.",
                sample_name,
                words.len(),
                words_per_entry
            )));
        }

        if words.len() / words_per_entry > u32::MAX as usize {
            return Err(JsValue::from_str(
                "Sample has too many split k-mers for u32 SNP distances.",
            ));
        }

        let progress_index = expected_index + 1;
        let on_stage = |stage: &str| {
            post_alignment_progress(
                self.run_id,
                stage,
                Some(progress_index),
                Some(self.sample_total),
                Some(&sample_name),
            );
        };

        logw(
            &format!(
                "Adding sample {} with {} split k-mers",
                sample_name,
                words.len() / words_per_entry
            ),
            Some("info"),
        );
        if self.k < 32 {
            self.alignment64
                .as_mut()
                .ok_or_else(|| JsValue::from_str("No 64-bit alignment is available."))?
                .add_sample_words(&sample_name, words, &on_stage)
                .map_err(|error| JsValue::from_str(&error))?;
        } else {
            self.alignment128
                .as_mut()
                .ok_or_else(|| JsValue::from_str("No 128-bit alignment is available."))?
                .add_sample_words(&sample_name, words, &on_stage)
                .map_err(|error| JsValue::from_str(&error))?;
        }
        self.file_names.push(sample_name);
        Ok(())
    }

    /// Build a tree from existing distances, including those of imported alignments.
    pub fn build_tree(&self) -> Result<String, JsValue> {
        self.tree_result(&|_| {}).map(|result| result.dump())
    }

    /// Finish this batch and rebuild its tree; alignment export remains lazy.
    pub fn finish_alignment(&mut self) -> Result<String, JsValue> {
        let can_export_alignment = if self.k < 32 {
            self.alignment64
                .as_ref()
                .is_some_and(SkaAlign::can_export_alignment)
        } else {
            self.alignment128
                .as_ref()
                .is_some_and(SkaAlign::can_export_alignment)
        };
        if !can_export_alignment {
            return Err(JsValue::from_str(
                "Batch finalisation requires an alignment built from split k-mers.",
            ));
        }
        if self.file_names.len() != self.sample_total {
            return Err(JsValue::from_str("All planned samples must be aggregated before finalisation."));
        }
        let on_stage = |stage: &str| {
            post_alignment_progress(self.run_id, stage, None, Some(self.sample_total), None);
        };
        let mut results = self.tree_result(&on_stage)?;
        results["alignmentAvailable"] = (self.file_names.len() >= 3).into();
        post_alignment_progress(
            self.run_id,
            "exporting",
            None,
            Some(self.sample_total),
            None,
        );
        Ok(results.dump())
    }

    /// Start a resumable export, splitting gzip members at the requested sample count.
    pub fn begin_alignment_export(&mut self, samples_per_part: u32) -> Result<(), JsValue> {
        if self.file_names.len() < 3 || self.file_names.len() != self.sample_total {
            return Err(JsValue::from_str(
                "At least three completed samples are required for alignment export.",
            ));
        }
        if samples_per_part == 0 {
            return Err(JsValue::from_str("Alignment export part size must be positive."));
        }
        logw("Obtaining alignment...", Some("info"));
        let result = if self.k < 32 {
            self.alignment64
                .as_mut()
                .ok_or_else(|| JsValue::from_str("No extracted alignment is available."))?
                .begin_alignment_export(samples_per_part as usize)
        } else {
            self.alignment128
                .as_mut()
                .ok_or_else(|| JsValue::from_str("No extracted alignment is available."))?
                .begin_alignment_export(samples_per_part as usize)
        };
        result.map_err(|error| JsValue::from_str(&error))?;
        post_alignment_progress(
            self.run_id,
            "obtaining-alignment",
            Some(0),
            Some(self.sample_total),
            None,
        );
        Ok(())
    }

    /// Return one compressed sample batch. Each returned byte vector is a gzip member.
    pub fn next_alignment_export_chunk(&mut self) -> Result<Vec<u8>, JsValue> {
        let run_id = self.run_id;
        let total = self.sample_total;
        let on_sample = |index| {
            post_alignment_progress(
                run_id,
                "obtaining-alignment",
                Some(index),
                Some(total),
                None,
            )
        };
        let (bytes, complete) = if self.k < 32 {
            let alignment = self
                .alignment64
                .as_mut()
                .ok_or_else(|| JsValue::from_str("No extracted alignment is available."))?;
            let bytes = alignment
                .next_alignment_export_chunk(&on_sample)
                .map_err(|error| JsValue::from_str(&error))?;
            (bytes, alignment.alignment_export_complete())
        } else {
            let alignment = self
                .alignment128
                .as_mut()
                .ok_or_else(|| JsValue::from_str("No extracted alignment is available."))?;
            let bytes = alignment
                .next_alignment_export_chunk(&on_sample)
                .map_err(|error| JsValue::from_str(&error))?;
            (bytes, alignment.alignment_export_complete())
        };
        if complete {
            logw(
                "Alignment chunks generated; awaiting output confirmation",
                Some("info"),
            );
        }
        Ok(bytes)
    }

    pub fn alignment_export_samples_processed(&self) -> usize {
        if self.k < 32 {
            self.alignment64
                .as_ref()
                .map_or(0, SkaAlign::alignment_export_samples_processed)
        } else {
            self.alignment128
                .as_ref()
                .map_or(0, SkaAlign::alignment_export_samples_processed)
        }
    }

    pub fn alignment_export_part_complete(&self) -> bool {
        if self.k < 32 {
            self.alignment64
                .as_ref()
                .is_some_and(SkaAlign::alignment_export_part_complete)
        } else {
            self.alignment128
                .as_ref()
                .is_some_and(SkaAlign::alignment_export_part_complete)
        }
    }

    pub fn alignment_export_complete(&self) -> bool {
        if self.k < 32 {
            self.alignment64
                .as_ref()
                .is_some_and(SkaAlign::alignment_export_complete)
        } else {
            self.alignment128
                .as_ref()
                .is_some_and(SkaAlign::alignment_export_complete)
        }
    }

    /// Release the sparse sites only after the final output was accepted by the caller.
    pub fn finish_alignment_export(&mut self) -> Result<(), JsValue> {
        let result = if self.k < 32 {
            self.alignment64
                .as_mut()
                .ok_or_else(|| JsValue::from_str("No extracted alignment is available."))?
                .finish_alignment_export()
        } else {
            self.alignment128
                .as_mut()
                .ok_or_else(|| JsValue::from_str("No extracted alignment is available."))?
                .finish_alignment_export()
        };
        result.map_err(|error| JsValue::from_str(&error))?;
        logw("Alignment export accepted; sparse sites released", Some("info"));
        Ok(())
    }

    /// Cancel the current output and reset the sparse cursors for a fresh export.
    pub fn cancel_alignment_export(&mut self) {
        if self.k < 32 {
            if let Some(alignment) = self.alignment64.as_mut() {
                alignment.cancel_alignment_export();
            }
        } else if let Some(alignment) = self.alignment128.as_mut() {
            alignment.cancel_alignment_export();
        }
    }

    /// Returns sample names (crate-internal, not exposed to JS).
    pub(crate) fn names(&self) -> &[String] {
        &self.file_names
    }

    /// Returns strict lower-triangle distances (crate-internal, not exposed to JS).
    pub(crate) fn flat_distances(&self) -> Result<&[u32], JsValue> {
        if let Some(alignment) = &self.alignment64 {
            Ok(alignment.get_flat_distances())
        } else if let Some(alignment) = &self.alignment128 {
            Ok(alignment.get_flat_distances())
        } else {
            Err(JsValue::from_str("No alignment distances are available."))
        }
    }

    /// Whether an imported pairwise distance was capped at u32::MAX.
    pub fn distances_capped(&self) -> bool {
        self.alignment64
            .as_ref()
            .is_some_and(SkaAlign::distances_capped)
            || self
                .alignment128
                .as_ref()
                .is_some_and(SkaAlign::distances_capped)
    }

    /// Returns graph JSON for the given SNP threshold.
    ///
    /// Format: `{"nodes":[{"id":"sample","cluster":1},...],
    ///           "links":[{"source":"s1","target":"s2","snp_distance":3.0},...]}`.
    pub fn get_graph_json(&self, threshold: f64) -> Result<String, JsValue> {
        let cutoff =
            checked_distance_cutoff(threshold).map_err(|error| JsValue::from_str(&error))?;
        let distances = self.flat_distances()?;
        if distances.is_empty() {
            return Ok("{\"nodes\":[],\"links\":[]}".to_string());
        }
        let (cluster_map, graph) = cluster_distances_flat(&self.file_names, distances, cutoff)
            .map_err(|error| JsValue::from_str(&error))?;

        let mut nodes = json::JsonValue::new_array();
        for name in &self.file_names {
            let mut node = json::JsonValue::new_object();
            node["id"] = name.clone().into();
            node["cluster"] = cluster_map[name].into();
            let _ = nodes.push(node);
        }

        let mut links = json::JsonValue::new_array();
        for edge in graph.edge_references() {
            let si = edge.source().index();
            let ti = edge.target().index();
            let (row, column) = if si < ti { (ti, si) } else { (si, ti) };
            let flat_idx = crate::wasm::alignment_matrix::strict_triangle_index(row, column)
                .ok_or_else(|| JsValue::from_str("Invalid strict-triangle distance index."))?;
            let dist = distances[flat_idx];
            let mut link = json::JsonValue::new_object();
            link["source"] = self.file_names[si].clone().into();
            link["target"] = self.file_names[ti].clone().into();
            link["snp_distance"] = dist.into();
            let _ = links.push(link);
        }

        let mut result = json::JsonValue::new_object();
        result["nodes"] = nodes;
        result["links"] = links;
        Ok(result.dump())
    }

    /// Returns pairwise SNP distances as a gzip-compressed square CSV matrix.
    pub fn get_distances_csv_gzip(&self) -> Result<Vec<u8>, JsValue> {
        let distances = self.flat_distances()?;
        distances_csv_gzip(&self.file_names, distances).map_err(|error| {
            JsValue::from_str(&format!("Could not gzip distance CSV output: {error}"))
        })
    }
}

#[cfg(test)]
mod alignment_export_tests {
    use super::{add_imported_mismatches, distances_csv_gzip};
    use flate2::read::GzDecoder;
    use std::io::Read;

    #[test]
    fn gzip_distance_csv_preserves_square_matrix_and_integer_values() {
        let names = vec![
            "sample-a".to_string(),
            "sample-b".to_string(),
            "sample-c".to_string(),
            "sample-d".to_string(),
        ];
        let compressed = distances_csv_gzip(&names, &[1, 2, 3, 4, 5, 6]).unwrap();
        let mut decoded = String::new();
        GzDecoder::new(compressed.as_slice())
            .read_to_string(&mut decoded)
            .unwrap();

        assert_eq!(
            decoded,
            ",sample-a,sample-b,sample-c,sample-d\n\
sample-a,0,1,2,4\n\
sample-b,1,0,3,5\n\
sample-c,2,3,0,6\n\
sample-d,4,5,6,0\n"
        );
    }

    #[test]
    fn gzip_distance_csv_is_empty_when_alignment_is_too_small() {
        assert!(
            distances_csv_gzip(&["one".to_string(), "two".to_string()], &[])
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn gzip_distance_csv_rejects_a_malformed_triangle() {
        let names = ["a", "b", "c", "d"].map(str::to_owned);
        assert!(distances_csv_gzip(&names, &[1, 2, 3]).is_err());
    }

    #[test]
    fn imported_distances_saturate_at_u32_max_and_report_it() {
        assert_eq!(add_imported_mismatches(u32::MAX - 1, 1), (u32::MAX, false));
        assert_eq!(add_imported_mismatches(u32::MAX, 1), (u32::MAX, true));
    }
}

#[cfg(test)]
mod fastq_pair_info_tests {
    use super::{
        fastq_files_are_mates, fastq_pair_info, prepare_alignment_samples, AlignmentSample,
    };

    #[test]
    fn pairs_only_when_the_read_mate_suffix_and_sample_match() {
        let mate_1 = fastq_pair_info("BARNARDS_Site_1_2017_ERR4920474_1.fastq.gz").unwrap();
        let mate_2 = fastq_pair_info("BARNARDS_Site_1_2017_ERR4920474_2.fastq.gz").unwrap();
        let different_sample =
            fastq_pair_info("BARNARDS_Site_7_2017_ERR4920592_2.fastq.gz").unwrap();

        assert_eq!(mate_1.0, mate_2.0);
        assert_eq!((mate_1.1, mate_2.1), (1, 2));
        assert_ne!(mate_1.0, different_sample.0);
        assert_eq!(mate_1.2, "BARNARDS_Site_1_2017_ERR4920474.fastq.gz");
        assert!(fastq_files_are_mates(
            "BARNARDS_Site_1_2017_ERR4920474_1.fastq.gz",
            "BARNARDS_Site_1_2017_ERR4920474_2.fastq.gz"
        ));
        assert!(!fastq_files_are_mates(
            "BARNARDS_Site_7_2017_ERR4920592_2.fastq.gz",
            "BARNARDS_Site_1_2017_ERR4920474_1.fastq.gz"
        ));
    }

    #[test]
    fn recognises_illumina_r_mates_with_chunk_suffixes() {
        let mate_1 = fastq_pair_info("sample_R1_001.fastq.gz").unwrap();
        let mate_2 = fastq_pair_info("sample_R2_001.fastq.gz").unwrap();

        assert_eq!(mate_1.0, mate_2.0);
        assert_eq!((mate_1.1, mate_2.1), (1, 2));
    }

    #[test]
    fn groups_mixed_inputs_and_preserves_alignment_sample_order() {
        let names = [
            "reads_R1.fastq.gz",
            "reference.fa",
            "BARNARDS_Site_1_2017_ERR4920474_1.fastq.gz",
            "BARNARDS_Site_1_2017_ERR4920474_2.fastq.gz",
            "BARNARDS_Site_7_2016_ERR4920603_1.fastq.gz",
            "BARNARDS_Site_7_2016_ERR4920603_2.fastq.gz",
            "BARNARDS_Site_7_2017_ERR4920592_1.fastq.gz",
            "BARNARDS_Site_7_2017_ERR4920592_2.fastq.gz",
        ]
        .map(str::to_owned);

        let samples = prepare_alignment_samples(&names);

        assert_eq!(samples.len(), 5);
        assert_eq!(samples[0].name, "reference.fa");
        assert_eq!(samples[1].name, "BARNARDS_Site_7_2017_ERR4920592.fastq.gz");
        assert_eq!(samples[1].first_index, 7);
        assert_eq!(samples[1].second_index, Some(6));
        assert_eq!(samples[2].name, "BARNARDS_Site_7_2016_ERR4920603.fastq.gz");
        assert_eq!(samples[2].first_index, 5);
        assert_eq!(samples[2].second_index, Some(4));
        assert_eq!(samples[3].name, "BARNARDS_Site_1_2017_ERR4920474.fastq.gz");
        assert_eq!(samples[3].first_index, 3);
        assert_eq!(samples[3].second_index, Some(2));
        assert_eq!(samples[4].name, "reads_R1.fastq.gz");
        assert_eq!(samples[4].first_index, 0);
    }

    #[test]
    fn two_unmatched_fastqs_are_two_samples_in_input_order() {
        let names = ["sample_a_1.fastq", "sample_b_2.fastq"].map(str::to_owned);

        let samples = prepare_alignment_samples(&names);

        assert_eq!(samples.len(), 2);
        assert_eq!(samples[0].first_index, 0);
        assert_eq!(samples[1].first_index, 1);
        assert!(samples.iter().all(|sample| sample.second_index.is_none()));
    }

    #[test]
    fn two_mate_fastqs_count_as_one_sample() {
        let names = ["sample_1.fastq.gz", "sample_2.fastq.gz"].map(str::to_owned);

        assert_eq!(
            prepare_alignment_samples(&names),
            vec![AlignmentSample {
                name: "sample.fastq.gz".to_owned(),
                first_index: 0,
                second_index: Some(1),
            }]
        );
    }
}
