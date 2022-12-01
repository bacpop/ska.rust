use clap::{Parser, Subcommand, ArgGroup};

fn valid_kmer(s: &str) -> Result<usize, String> {
    let k: usize = s
        .parse()
        .map_err(|_| format!("`{}` isn't a valid k-mer", s))?;
    if k < 5 || k > 31 || k % 2 == 0 {
        Err(format!("K-mer must an odd number between 5 and 31 (inclusive)"))
    } else {
        Ok(k)
    }
}

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
#[command(propagate_version = true)]
pub struct Args {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
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
        #[arg(short, value_parser = valid_kmer, default_value_t = 17)]
        k: usize,

        /// Ignore reverse complement (all contigs are oriented along same strand).
        #[arg(long, default_value_t = false)]
        single_strand: bool,

        /// Number of CPU threads
        #[arg(long, default_value_t = 1)]
        threads: usize,
     },
    // ska merge <files.skf>... -o <prefix>
    /// Combine multiple split k-mer files
    Merge { name: Option<String> },
    // ska delete <file.skf> -l <name_list>
    /// Remove samples from a split k-mer file
    Delete { name: Option<String> },
    // ska weed <file.skf> --seqs <list.fa>
    /// Remove k-mers from a split k-mer file
    Weed { name: Option<String> },
    // ska align <input_files>... <-o output.aln> [--min <freq>] [--const-sites] [--cpus <cpus>]
    /// Write an unordered alignment
    Align { name: Option<String> },
    // ska map <ref.fa> <input_files>... <-o output.aln> [--min <freq>] [--const-sites] [--cpus <cpus>]
    /// Write an ordered alignment using a reference sequence
    Map { name: Option<String> },
    // ska nk file.skf <--full-info>
    /// Get the number of k-mers in a split k-mer file, and other information
    Nk { name: Option<String> },

}

/*
--seqs <list.fa>  A fasta file with sequences to remove.
--min <freq>      Minimum frequency present to output to alignment [default: 0.9].
--const-sites     Write invariant sites.
--full-info       Write out split-kmers, and middle base matrix.

--cpus <cpus>     Use parallel processing, with this many cores.
*/

pub fn cli_args() -> Args {
    let args = Args::parse();
    return args;
}