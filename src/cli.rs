use docopt::Docopt;
use serde::Deserialize;

const USAGE: &'static str = "
Usage:
    ska build <seq_files>... -o <prefix> [-k <k-size>] [--single-strand] [--compress] [--cpus <cpus>]
    ska build -l <file_list> -o <prefix> [-k <k-size>] [--single-strand] [--compress] [--cpus <cpus>]

    ska merge <files.skf>... -o <prefix>
    ska delete <file.skf> -l <name_list>
    ska weed <file.skf> --seqs <list.fa>

    ska align <seq_files>... <-o output.aln> [--min <freq>] [--const-sites] [--cpus <cpus>]
    ska map <ref.fa> <seq_files>... <-o output.aln> [--min <freq>] [--const-sites] [--cpus <cpus>]

    ska nk file.skf <--full-info>

    ska (-h | --help)
    ska (--version)

Options:
    -h --help         Show this help.
    --version         Show version.

    -o <prefix>       Output prefix.
    -k <k-size>       K-mer size [default: 17].
    --single-strand   Ignore reverse complement (all contigs are oriented along same strand).
    --compress        Compress the output skf file.
    -l <file_list>    A list of input files, one per line, tab separated name and file.
    --seqs <list.fa>  A fasta file with sequences to remove.
    --min <freq>      Minimum frequency present to output to alignment [default: 0.9].
    --const-sites     Write invariant sites.
    --full-info       Write out split-kmers, and middle base matrix.

    --cpus <cpus>     Use parallel processing, with this many cores.
";

#[derive(Deserialize)]
struct Args {
    flag_build: bool,
    flag_merge: bool,
    flag_delete: bool,
    flag_weed: bool,
    flag_align: bool,
    flag_map: bool,
    flag_nk: bool,
    flag_h: bool,
    flag_help: bool,
    flag_version: bool,

    arg_prefix: String,
    arg_k: usize,
    arg_single_strand: bool,
    arg_compress: bool,
    arg_seq_files: Vec<String>,
    arg_file_list: String,
    arg_seqs: String,
    arg_min: f64,
    arg_const_sites: bool,
    arg_full_info: bool,
    arg_cpus: usize,
}

fn cli_args() -> Args {
    let argv = std::env::args();
    let args: Args = Docopt::new(USAGE)
        .and_then(|d| d.argv(argv().into_iter()).deserialize())
        .unwrap_or_else(|e| e.exit());

    if args.flag_version {
        // TODO work out how to do this properly
    }

    return args;
}