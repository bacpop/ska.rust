//! Command line interface, built using [`crate::clap` with `Derive`](https://docs.rs/clap/latest/clap/_derive/_tutorial/index.html)
use std::fmt;

use clap::{ArgGroup, Parser, Subcommand, ValueEnum};

use super::QualFilter;

/// Default split k-mer size
pub const DEFAULT_KMER: usize = 17;
/// Default single strand (which is equivalent to !rc)
pub const DEFAULT_STRAND: bool = false;
/// Default minimum k-mer count for FASTQ files
pub const DEFAULT_MINCOUNT: u16 = 10;
/// Default minimum base quality (PHRED score) for FASTQ files
pub const DEFAULT_MINQUAL: u8 = 20;
/// Default quality filtering criteria
pub const DEFAULT_QUALFILTER: QualFilter = QualFilter::Middle;

#[doc(hidden)]
fn valid_kmer(s: &str) -> Result<usize, String> {
    let k: usize = s
        .parse()
        .map_err(|_| format!("`{s}` isn't a valid k-mer"))?;
    if !(5..=63).contains(&k) || k % 2 == 0 {
        Err("K-mer must an odd number between 5 and 63 (inclusive)".to_string())
    } else {
        Ok(k)
    }
}

#[doc(hidden)]
fn zero_to_one(s: &str) -> Result<f64, String> {
    let f: f64 = s
        .parse()
        .map_err(|_| format!("`{s}` isn't a valid frequency"))?;
    if !(0.0..=1.0).contains(&f) {
        Err("Frequency must be between 0 and 1 (inclusive)".to_string())
    } else {
        Ok(f)
    }
}

#[doc(hidden)]
fn valid_cpus(s: &str) -> Result<usize, String> {
    let threads: usize = s
        .parse()
        .map_err(|_| format!("`{s}` isn't a valid number of cores"))?;
    if threads < 1 {
        Err("Threads must be one or higher".to_string())
    } else {
        Ok(threads)
    }
}

/// Prints a warning if more threads than available have been requested
pub fn check_threads(threads: usize) {
    let max_threads = num_cpus::get();
    if threads > max_threads {
        log::warn!("{threads} threads is greater than available cores {max_threads}");
    }
}

/// Possible output file types
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
pub enum FileType {
    /// Variant call format
    Vcf,
    /// FASTA alignment
    Aln,
}

/// Possible variant filters for align
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
pub enum FilterType {
    /// Output all variants
    NoFilter,
    /// Filter constant bases
    NoConst,
    /// Filter constant bases, and any variants where the only variation is an ambiguous base
    NoAmbigOrConst,
}

/// As text, for use in logging messages
impl fmt::Display for FilterType {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            Self::NoFilter => write!(f, "No filtering"),
            Self::NoConst => write!(f, "No constant sites"),
            Self::NoAmbigOrConst => write!(f, "No constant sites or ambiguous bases"),
        }
    }
}

/// Options that apply to all subcommands
#[derive(Parser)]
#[command(author, version, about, long_about = None)]
#[command(propagate_version = true)]
pub struct Args {
    #[doc(hidden)]
    #[command(subcommand)]
    pub command: Commands,

    /// Show progress messages
    #[arg(short, long, global = true)]
    pub verbose: bool,
}

/// Subcommands and their specific options
#[derive(Subcommand)]
pub enum Commands {
    #[command(group(
        ArgGroup::new("input")
            .required(true)
            .args(["seq_files", "file_list"]),
    ))]
    /// Create a split-kmer file from input sequences
    Build {
        /// List of input FASTA files
        #[arg(group = "input")]
        seq_files: Option<Vec<String>>,

        /// File listing input files (tab separated name, sequences)
        #[arg(short, group = "input")]
        file_list: Option<String>,

        /// Output prefix
        #[arg(short)]
        output: String,

        /// K-mer size
        #[arg(short, value_parser = valid_kmer, default_value_t = DEFAULT_KMER)]
        k: usize,

        /// Ignore reverse complement (all contigs are oriented along same strand)
        #[arg(long, default_value_t = DEFAULT_STRAND)]
        single_strand: bool,

        /// Minimum k-mer count (with reads)
        #[arg(long, default_value_t = DEFAULT_MINCOUNT)]
        min_count: u16,

        /// Minimum k-mer count (with reads)
        #[arg(long, default_value_t = DEFAULT_MINQUAL)]
        min_qual: u8,

        /// Quality filtering criteria (with reads)
        #[arg(long, value_enum, default_value_t = QualFilter::Middle)]
        qual_filter: QualFilter,

        /// Number of CPU threads
        #[arg(long, value_parser = valid_cpus, default_value_t = 1)]
        threads: usize,
    },
    /// Write an unordered alignment
    Align {
        /// A .skf file, or list of .fasta files
        #[arg(required = true)]
        input: Vec<String>,

        /// Output filename (omit to output to stdout)
        #[arg(short)]
        output: Option<String>,

        /// Minimum fraction of samples a k-mer has to appear in
        #[arg(short, long, value_parser = zero_to_one, default_value_t = 0.9)]
        min_freq: f64,

        /// Filter for constant middle base sites
        #[arg(long, value_enum, default_value_t = FilterType::NoConst)]
        filter: FilterType,

        /// Number of CPU threads
        #[arg(long, value_parser = valid_cpus, default_value_t = 1)]
        threads: usize,
    },
    /// Write an ordered alignment using a reference sequence
    Map {
        /// Reference FASTA file to map to
        reference: String,

        /// A .skf file, or list of .fasta files
        input: Vec<String>,

        /// Output filename (omit to output to stdout)
        #[arg(short)]
        output: Option<String>,

        /// Format of output file
        #[arg(short, long, value_enum, default_value_t = FileType::Aln)]
        format: FileType,

        /// Number of CPU threads
        #[arg(long, value_parser = valid_cpus, default_value_t = 1)]
        threads: usize,
    },
    /// Calculate SNP distances and single-linkage clusters
    Distance {
        /// Split-kmer (.skf) file to operate on
        skf_file: String,

        /// Output filename (omit to output to stdout)
        #[arg(short)]
        output: Option<String>,

        /// Number of CPU threads
        #[arg(long, value_parser = valid_cpus, default_value_t = 1)]
        threads: usize,
    },
    /// Combine multiple split k-mer files
    Merge {
        /// List of input split-kmer (.skf) files
        skf_files: Vec<String>,

        /// Output prefix
        #[arg(short)]
        output: String,
    },
    #[command(group(
        ArgGroup::new("input")
            .required(true)
            .args(["names", "file_list"]),
    ))]
    /// Remove samples from a split k-mer file
    Delete {
        /// Split-kmer (.skf) file to operate on
        skf_file: String,

        /// File listing sample names to remove
        #[arg(short, group = "input")]
        file_list: Option<String>,

        /// List of sample names to remove
        #[arg(group = "input")]
        names: Option<Vec<String>>,
    },
    /// Remove k-mers from a split k-mer file
    Weed {
        /// Split-kmer (.skf) file to operate on
        skf_file: String,

        /// A FASTA file containing sequences to remove
        weed_file: Option<String>,

        /// Remove k-mers not in the weed_file
        #[arg(long, default_value_t = false)]
        reverse: bool,

        /// Minimum fraction of samples a k-mer has to appear in
        #[arg(short, long, value_parser = zero_to_one, default_value_t = 0.0)]
        min_freq: f64,

        /// Filter for constant middle base sites
        #[arg(long, value_enum, default_value_t = FilterType::NoFilter)]
        filter: FilterType,
    },
    /// Get the number of k-mers in a split k-mer file, and other information
    Nk {
        /// Split-kmer (.skf) file to operate on
        skf_file: String,

        /// Also write out split-kmers, and middle base matrix
        #[arg(long, default_value_t = false)]
        full_info: bool,
    },
    /// Estimate a coverage cutoff using a k-mer count profile (FASTQ only)
    Cov {
        /// FASTQ file (or .fastq.gz) with forward reads
        fastq_fwd: String,

        /// FASTQ file (or .fastq.gz) with reverse reads
        fastq_rev: String,

        /// K-mer size
        #[arg(short, value_parser = valid_kmer, default_value_t = DEFAULT_KMER)]
        k: usize,

        /// Ignore reverse complement (all reads are oriented along same strand)
        #[arg(long, default_value_t = DEFAULT_STRAND)]
        single_strand: bool,
    },
}

/// Function to parse command line args into [`Args`] struct
pub fn cli_args() -> Args {
    Args::parse()
}
