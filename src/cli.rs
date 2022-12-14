use clap::{ArgGroup, Parser, Subcommand, ValueEnum};

extern crate num_cpus;

pub const DEFAULT_KMER: usize = 17;
pub const DEFAULT_STRAND: bool = false;
pub const DEFAULT_MINCOUNT: u16 = 10;
pub const DEFAULT_MINQUAL: u8 = 20;

fn valid_kmer(s: &str) -> Result<usize, String> {
    let k: usize = s
        .parse()
        .map_err(|_| format!("`{}` isn't a valid k-mer", s))?;
    if k < 5 || k > 31 || k % 2 == 0 {
        Err(format!(
            "K-mer must an odd number between 5 and 31 (inclusive)"
        ))
    } else {
        Ok(k)
    }
}

fn zero_to_one(s: &str) -> Result<f64, String> {
    let f: f64 = s
        .parse()
        .map_err(|_| format!("`{}` isn't a valid frequency", s))?;
    if f < 0.0 || f > 1.0 {
        Err(format!("Frequency must be between 0 and 1 (inclusive)"))
    } else {
        Ok(f)
    }
}

fn valid_cpus(s: &str) -> Result<usize, String> {
    let threads: usize = s
        .parse()
        .map_err(|_| format!("`{}` isn't a valid number of cores", s))?;
    let max_threads = num_cpus::get();
    if threads < 1 || threads > max_threads {
        Err(format!("Threads must be between 1 and {}", max_threads))
    } else {
        Ok(threads)
    }
}

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
pub enum FileType {
    /// Variant call format
    Vcf,
    /// FASTA alignment
    Aln,
}

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
#[command(propagate_version = true)]
pub struct Args {
    #[command(subcommand)]
    pub command: Commands,

    /// Show progress messages
    #[arg(short, long, global = true)]
    pub verbose: bool,
}

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

        /// Number of CPU threads
        #[arg(long, value_parser = valid_cpus, default_value_t = 1)]
        threads: usize,
    },
    /// Write an unordered alignment
    Align {
        /// A .skf file, or list of .fasta files
        #[arg(required = true)]
        input: Vec<String>,

        /// Output prefix (omit to output to stdout)
        #[arg(short)]
        output: Option<String>,

        /// Minimum fraction of samples a k-mer has to appear in
        #[arg(short, long, value_parser = zero_to_one, default_value_t = 0.9)]
        min_freq: f64,

        /// Output constant middle base sites
        #[arg(long, default_value_t = false)]
        const_sites: bool,

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

        /// Output prefix
        #[arg(short)]
        output: Option<String>,

        #[arg(short, long, value_enum, default_value_t = FileType::Aln)]
        format: FileType,

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

        /// Minimum fraction of samples a k-mer has to appear in
        #[arg(short, long, value_parser = zero_to_one, default_value_t = 0.0)]
        min_freq: f64,

        /// Remove constant middle base sites
        #[arg(long, default_value_t = false)]
        remove_const_sites: bool,
    },
    /// Get the number of k-mers in a split k-mer file, and other information
    Nk {
        /// Split-kmer (.skf) file to operate on
        skf_file: String,

        /// Also write out split-kmers, and middle base matrix
        #[arg(long, default_value_t = false)]
        full_info: bool,
    },
}

pub fn cli_args() -> Args {
    let args = Args::parse();
    return args;
}
