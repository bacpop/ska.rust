//! Split k-mer analysis (version 2)  uses exact matching of split k-mer sequences to align closely related
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
//! K-mer size must be odd, greater than 5, and a maximum of 31. Smaller k-mers
//! are more sensitive and can find closer positions, but are less specific so
//! may lead to more repeated split k-mers and ambiguous bases.
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
//! To use FASTQ files, or specify sample names, or more easily input a larger number of input files,
//! you can provide a tab separated file list via `-f` instead of listing files. For example
//! a file called `input_sequence.txt` to describe `sample_1` (an assembly) and `sample_2` (paired reads):
//! ```text
//! sample_1    assemblies/seq1.fa
//! sample_2    reads/seq2_1.fastq.gz   reads/seq2_2.fastq.gz
//! ```
//! You'd run:
//! ```bash
//! ska build -o seqs -f input_sequence.txt --min-count 20 --min-qual 30
//! ```
//! (here also demonstrating changing the error filtering criteria with the FASTQ files.)
//!
//! FASTQ files must be paired end. If you'd like to request more flexibility in
//! this regard please contact us.
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
//! ska align --min-freq 1 --filter NoFilter -o seqs seqs.skf
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
//! Create an alignment from a `.skf` file or sequence files by mapping its
//! split k-mers to split k-mers of a linear reference sequence. This produces pseudoalignment
//! with the same sites/columns as the reference sequence. Sites not mapped will
//! be output as missing '-'.
//!
//! Pass the FASTA file of the reference as the first argument, and the skf file as
//! the second argument:
//! ```bash
//! ska map ref.fa seqs.skf -o ref_mapped.aln
//! ```
//!
//! You can also get a VCF as output, which has rows as variants, and only has the
//! variable sites (but will include unmapped bases as missing).
//! An example command, also demonstrating that everything can be done from input sequences in a single command:
//! ```bash
//! ska map ref.fa seq1.fa seq2.fa -f vcf --threads 2 | bgzip -c - > seqs.vcf.gz
//! ```
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
//! file, or an error will be thrown. After removing samples, the input `.skf` file will be overwritten.
//! Any split k-mers which have no observed missing bases after removal will also be deleted.
//!
//! ```bash
//! ska delete all_samples.skf sample_1 sample_3
//! ```
//! If you wish to delete many samples, you can use `-f` to provide their names
//! by line via an input file.
//!
//! ## ska weed
//!
//! Remove (weed) split k-mers from an `.skf` file. The split k-mers to be removed
//! are generated from a FASTA file, which may for example contain known transposons or
//! other repeat sequences. As with `delete`, the input `.skf` file is overwritten.
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
//! ska weed --filter NoAmbigOrConst --min-freq 0.9 all_samples.skf
//! ```
//!
//! ## ska nk
//! Return information on the *n*umber of *k*-mers in an `.skf` file. This will
//! print on STDOUT:
//! - The k-mer size.
//! - Whether reverse complements were used.
//! - The number of split k-mers.
//! - The number of samples.
//! - The sample names.
//!
//! ```bash
//! ska nk all_samples.skf
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
//! # API usage
//!
//! See the submodule documentation linked below.
//!
//! To use `ska build` in other rust code:
//! ```rust
//! use ska::merge_ska_dict::{InputFastx, build_and_merge};
//! use ska::merge_ska_array::MergeSkaArray;
//!
//! // Build, merge
//! let input_files: [InputFastx; 2] = [("test1".to_string(),
//!                                      "tests/test_files_in/test_1.fa".to_string(),
//!                                      None),
//!                                     ("test2".to_string(),
//!                                      "tests/test_files_in/test_2.fa".to_string(),
//!                                      None)];
//! let rc = true;
//! let k = 17;
//! let min_count = 1;
//! let min_qual = 0;
//! let threads = 2;
//! let merged_dict =
//!     build_and_merge(&input_files, k, rc, min_count, min_qual, threads);
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
//! let mut ska_array = load_array(&input, threads);
//!
//! // Apply filters
//! let min_count = 2;
//! let update_kmers = false;
//! let filter = FilterType::NoConst;
//! ska_array.filter(min_count, &filter, update_kmers);
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
use std::cmp::min;
use std::time::Instant;

pub mod merge_ska_dict;
pub mod ska_dict;
use crate::merge_ska_dict::build_and_merge;

pub mod ska_ref;
use crate::ska_ref::RefSka;
pub mod merge_ska_array;
use crate::merge_ska_array::MergeSkaArray;

pub mod generic_modes;
use crate::generic_modes::*;

pub mod cli;
use crate::cli::*;

pub mod io_utils;
use crate::io_utils::*;

#[doc(hidden)]
pub fn main() {
    let args = cli_args();
    if args.verbose {
        simple_logger::init_with_level(log::Level::Info).unwrap();
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
            threads,
        } => {
            // Read input
            let input_files = get_input_list(file_list, seq_files);

            // Build, merge
            let rc = !*single_strand;
            if *k <= 31 {
                log::info!("k={}: using 64-bit representation", *k);
                type IntSize = u64;
                let merged_dict = build_and_merge::<IntSize>(
                    &input_files,
                    *k,
                    rc,
                    *min_count,
                    *min_qual,
                    *threads,
                );

                // Save
                log::info!("Converting to array representation and saving");
                let ska_array = MergeSkaArray::new(&merged_dict);
                ska_array
                    .save(format!("{output}.skf").as_str())
                    .expect("Failed to save output file");
            } else {
                log::info!("k={}: using 128-bit representation", *k);
                type IntSize = u128;
                let merged_dict = build_and_merge::<IntSize>(
                    &input_files,
                    *k,
                    rc,
                    *min_count,
                    *min_qual,
                    *threads,
                );

                // Save
                log::info!("Converting to array representation and saving");
                let ska_array = MergeSkaArray::new(&merged_dict);
                ska_array
                    .save(format!("{output}.skf").as_str())
                    .expect("Failed to save output file");
            }
        }
        Commands::Align {
            input,
            output,
            min_freq,
            filter,
            threads,
        } => {
            if let Ok(mut ska_array) = load_array::<u64>(input, *threads) {
                // In debug mode (cannot be set from CLI, give details)
                log::debug!("{ska_array}");
                align(&mut ska_array, output, filter, *min_freq);
            } else if let Ok(mut ska_array) = load_array::<u128>(input, *threads) {
                // In debug mode (cannot be set from CLI, give details)
                log::debug!("{ska_array}");
                align(&mut ska_array, output, filter, *min_freq);
            } else {
                panic!("Could not read input file(s): {:?}", input);
            }
        }
        Commands::Map {
            reference,
            input,
            output,
            format,
            threads,
        } => {
            log::info!("Loading skf as dictionary");
            if let Ok(mut ska_array) = load_array::<u64>(input, *threads) {
                log::info!(
                    "Making skf of reference k={} rc={}",
                    ska_array.kmer_len(),
                    ska_array.rc()
                );
                let mut ska_ref =
                    RefSka::<u64>::new(ska_array.kmer_len(), reference, ska_array.rc());
                map(&mut ska_array, &mut ska_ref, output, format, *threads);
            } else if let Ok(mut ska_array) = load_array::<u128>(input, *threads) {
                log::info!(
                    "Making skf of reference k={} rc={}",
                    ska_array.kmer_len(),
                    ska_array.rc()
                );
                let mut ska_ref =
                    RefSka::<u128>::new(ska_array.kmer_len(), reference, ska_array.rc());
                map(&mut ska_array, &mut ska_ref, output, format, *threads);
            } else {
                panic!("Could not read input file(s): {:?}", input);
            }
        }
        Commands::Merge { skf_files, output } => {
            if skf_files.len() < 2 {
                panic!("Need at least two files to merge");
            }

            log::info!("Loading first alignment");
            if let Ok(mut first_array) = MergeSkaArray::<u64>::load(&skf_files[0]) {
                merge(&mut first_array, &skf_files[1..], output);
            } else if let Ok(mut first_array) = MergeSkaArray::<u128>::load(&skf_files[0]) {
                merge(&mut first_array, &skf_files[1..], output);
            } else {
                panic!("Could not read input file: {}", skf_files[0]);
            }
        }
        Commands::Delete {
            skf_file,
            file_list,
            names,
        } => {
            let input_files = get_input_list(file_list, names);
            let input_names: Vec<&str> = input_files.iter().map(|t| &*t.0).collect();
            log::info!("Loading skf file");
            if let Ok(mut ska_array) = MergeSkaArray::<u64>::load(&skf_file) {
                delete(&mut ska_array, &input_names, skf_file);
            } else if let Ok(mut ska_array) = MergeSkaArray::<u128>::load(&skf_file) {
                delete(&mut ska_array, &input_names, skf_file);
            } else {
                panic!("Could not read input file: {skf_file}");
            }
        }
        Commands::Weed {
            skf_file,
            weed_file,
            min_freq,
            filter,
        } => {
            log::info!("Loading skf file");
            if let Ok(mut ska_array) = MergeSkaArray::<u64>::load(&skf_file) {
                weed(&mut ska_array, weed_file, *min_freq, filter, &skf_file);
            } else if let Ok(mut ska_array) = MergeSkaArray::<u128>::load(&skf_file) {
                weed(&mut ska_array, weed_file, *min_freq, filter, &skf_file);
            } else {
                panic!("Could not read input file: {skf_file}");
            }
        }
        Commands::Nk {
            skf_file,
            full_info,
        } => {
            if let Ok(ska_array) = MergeSkaArray::<u64>::load(&skf_file) {
                println!("{ska_array}");
                if *full_info {
                    log::info!("Printing full info");
                    println!("{ska_array:?}");
                }
            } else if let Ok(ska_array) = MergeSkaArray::<u128>::load(&skf_file) {
                println!("{ska_array}");
                if *full_info {
                    log::info!("Printing full info");
                    println!("{ska_array:?}");
                }
            } else {
                panic!("Could not read input file: {skf_file}");
            }
        }
    }
    let end = Instant::now();

    eprintln!("SKA done in {}s", end.duration_since(start).as_secs());
    eprintln!("⬛⬜⬛⬜⬛⬜⬛⬜⬛⬜");
    eprintln!("⬜⬛⬜⬛⬜⬛⬜⬛⬜⬛");
    eprintln!("⬛⬜⬛⬜⬛⬜⬛⬜⬛⬜");
    eprintln!("⬜⬛⬜⬛⬜⬛⬜⬛⬜⬛");
    log::info!("Complete");
}
