//! Common helper functions for parsing file input, loading, and setting output
//!
//! The functions are used by a few different subcommands to set correct
//! args to build structs, given the command line input

use std::error::Error;
use std::fs::File;
use std::io::{stdout, BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use regex::Regex;

use super::QualOpts;
use crate::merge_ska_array::MergeSkaArray;
use crate::merge_ska_dict::{build_and_merge, InputFastx};
use crate::ska_dict::bit_encoding::UInt;
use crate::CoverageHistogram;
use crate::cli::check_threads;

use crate::cli::{
    DEFAULT_KMER, DEFAULT_MINCOUNT, DEFAULT_MINQUAL, DEFAULT_PROPORTION_READS, DEFAULT_QUALFILTER,
    DEFAULT_STRAND,
};
use crate::ValidMinKmer;

/// Given a list of input files, parses them into triples of name, filename and
/// [`None`] to be used with [SkaDict](`crate::ska_dict::SkaDict::new()`).
///
/// To form the name, common file extensions are removed, and the basepath is
/// used. If this cannot be parsed then the full filename is used
pub fn read_input_fastas(seq_files: &[String]) -> Vec<InputFastx> {
    let mut input_files = Vec::new();
    // matches the file name (no extension) in a full path
    let re_path = Regex::new(r"^.+/(.+)\.(?i:fa|fasta|fastq|fastq\.gz)$").unwrap();
    // matches the file name (no extension) with no path
    let re_name = Regex::new(r"^(.+)\.(?i:fa|fasta|fastq|fastq\.gz)$").unwrap();
    for file in seq_files {
        let caps = re_path.captures(file).or(re_name.captures(file));
        let name = match caps {
            Some(capture) => capture[1].to_string(),
            None => file.to_string(),
        };
        input_files.push((name, file.to_string(), None));
    }
    input_files
}

/// Given a list of files via the CLI, loads or creates a
/// [MergeSkaArray](`crate::merge_ska_array::MergeSkaArray`).
///
/// This is to support functions like `ska align` and `ska map` working
/// from either a .skf file previously produced by `ska build`, or from a file
/// list so everything can be done in a single command.
///
/// If a single input file is given, then try and load it as an .skf file
///
/// If multiple files are given, they are assumed to be FASTA and
/// [build_and_merge](`crate::merge_ska_dict::build_and_merge`) is used to
/// create a merged dictionary, which is then converted to a merged array.
pub fn load_array<IntT: for<'a> UInt<'a>>(
    input: &[String],
    threads: usize,
) -> Result<MergeSkaArray<IntT>, Box<dyn Error>> {
    // Obtain a merged ska array
    if input.len() == 1 {
        log::info!(
            "Single file as input, trying to load as skf {}-bits",
            IntT::n_bits()
        );
        if threads > 1 {
            log::warn!("--threads only used if building skf, setting to 1");
            check_threads(1);
        }
        MergeSkaArray::load(input[0].as_str())
    } else {
        log::info!("Multiple files as input, running ska build with default settings");
        let input_files = read_input_fastas(input);
        let default_qual = QualOpts {
            min_count: DEFAULT_MINCOUNT,
            min_qual: DEFAULT_MINQUAL,
            qual_filter: DEFAULT_QUALFILTER,
        };
        let merged_dict = build_and_merge(
            &input_files,
            DEFAULT_KMER,
            !DEFAULT_STRAND,
            &default_qual,
            threads,
            DEFAULT_PROPORTION_READS,
        );
        Ok(MergeSkaArray::new(&merged_dict))
    }
}

/// Set a buffered stream to write to.
///
/// Either a file (if [`Some`]) or stdout otherwise (if [`None`]).
pub fn set_ostream(oprefix: &Option<String>) -> BufWriter<Box<dyn Write>> {
    let out_writer = match oprefix {
        Some(prefix) => {
            let path = Path::new(prefix);
            Box::new(File::create(path).unwrap()) as Box<dyn Write>
        }
        None => Box::new(stdout()) as Box<dyn Write>,
    };
    BufWriter::new(out_writer)
}

/// Obtain a list of input files and names from command line input.
///
/// If `file_list` is provided, read each line as `name\tseq1\tseq2`, where
/// `seq2` is optional, and if present the reverse fastqs. Otherwise, treat
/// as fasta.
///
/// If `seq_files` are provided use [`read_input_fastas`].
pub fn get_input_list(
    file_list: &Option<String>,
    seq_files: &Option<Vec<String>>,
) -> Vec<InputFastx> {
    // Read input
    match file_list {
        Some(files) => {
            let mut input_files: Vec<InputFastx> = Vec::new();
            let f = File::open(files).expect("Unable to open file_list");
            let f = BufReader::new(f);
            for line in f.lines() {
                let line = line.expect("Unable to read line in file_list");
                let fields: Vec<&str> = line.split_whitespace().collect();
                // Should be 2 entries for fasta, 3 for fastq
                let second_file = match fields.len() {
                    0..=1 => {
                        panic!("Unable to parse line in file_list")
                    }
                    2 => None,
                    3 => Some(fields[2].to_string()),
                    _ => {
                        panic!("Unable to parse line in file_list")
                    }
                };
                input_files.push((fields[0].to_string(), fields[1].to_string(), second_file));
            }
            input_files
        }
        None => read_input_fastas(seq_files.as_ref().unwrap()),
    }
}

/// Checks if any input files are fastq
pub fn any_fastq(files: &[InputFastx]) -> bool {
    files.iter().any(|file| file.2.is_some())
}

/// Counts number of fastq files
pub fn count_fastq(files: &[InputFastx]) -> usize {
    files.iter().filter(|file| file.2.is_some()).count()
}

/// Collects filepaths to the first 2 fastq files into a tuple
pub fn get_2_fastq_path(files: &[InputFastx]) -> (String, String) {
    let out: Vec<String> = files
        .iter()
        .filter(|file| file.2.is_some())
        .take(2)
        .map(|x| x.1.clone())
        .collect();

    if out.len().gt(&1) {
        (out[0].clone(), out[1].clone())
    } else {
        panic!("Trying to get 2 fastq files from a vector with <2 elements");
    }
}

/// Calculates minimum kmer cutoff depending on user provided argument
pub fn kmer_min_cutoff<IntT: for<'a> UInt<'a>>(
    v: &Option<ValidMinKmer>,
    files: &[InputFastx],
    k: &usize,
    rc: bool,
    verbose: bool,
) -> u16 {
    // Minimum kmer cutoff logic
    if v.is_none() {
        log::info!(
            "Using user-provided minimum kmer value of {}",
            DEFAULT_MINCOUNT
        );
        DEFAULT_MINCOUNT
    } else {
        match v.unwrap() {
            // User-provided value (already checked by cli parser)
            ValidMinKmer::Val(x) => {
                log::info!("Using user-provided minimum kmer value of {}", x);
                x
            }
            // auto-calculate & there are enough fastq files
            ValidMinKmer::Auto if count_fastq(files).ge(&2) => {
                let (fastq_fwd, fastq_rev) = get_2_fastq_path(files);
                let mut cov =
                    CoverageHistogram::<IntT>::new(&fastq_fwd, &fastq_rev, *k, rc, verbose);
                let out = cov.fit_histogram().expect("Couldn't fit coverage model") as u16;
                cov.plot_hist();
                log::info!("Using inferred minimum kmer value of {}", out);
                out
            }
            // Not enough fastq files, use default and warn user
            ValidMinKmer::Auto => {
                log::info!(
                    "Not enough fastq files to fit mixture model, using default kmer count of 5"
                );
                DEFAULT_MINCOUNT
            }
        }
    }
}
