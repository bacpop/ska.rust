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
use crate::ska_dict::bit_encoding::{decode_kmer, UInt};
use crate::skalo_utils::{encode_kmer, rev_compl};
use hashbrown::HashMap;

use crate::cli::{
    DEFAULT_KMER, DEFAULT_MINCOUNT, DEFAULT_MINQUAL, DEFAULT_PROPORTION_READS, DEFAULT_QUALFILTER,
    DEFAULT_STRAND,
};

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
    let mut fastq = false;
    for file in files {
        if file.2.is_some() {
            fastq = true;
            break;
        }
    }
    fastq
}

/// Reads .skf files for skalo
pub fn read_skalo_input(
    input_file: &str,
) -> (
    usize,
    Vec<String>,
    HashMap<u128, HashMap<u128, u32>>,
    HashMap<u32, String>,
) {
    println!(" # read file '{}'", input_file);

    // read the skf file and load split-kmers (ska_array), kmer length and sample names
    let ska_array = load_array::<u128>(&[input_file.to_string()], 1)
        .expect("\nerror: coud not read the skf file\n\n");
    let sample_names = ska_array.names().to_vec();
    let len_kmer = ska_array.kmer_len();
    let mask = (1 << (len_kmer * 2)) - 1;

    println!("     . {}-mers", len_kmer);
    println!("     . {} samples", sample_names.len());

    println!(" # build colored de Bruijn graph");

    // build De Bruijn graph
    let degenerate_code: HashMap<u8, Vec<char>> = [
        (b'A', vec!['A']),
        (b'T', vec!['T']),
        (b'G', vec!['G']),
        (b'C', vec!['C']),
        (b'M', vec!['A', 'C']),
        (b'S', vec!['C', 'G']),
        (b'W', vec!['A', 'T']),
        (b'R', vec!['A', 'G']),
        (b'Y', vec!['C', 'T']),
        (b'K', vec!['G', 'T']),
        (b'B', vec!['C', 'G', 'T']),
        (b'D', vec!['A', 'G', 'T']),
        (b'H', vec!['A', 'C', 'T']),
        (b'V', vec!['A', 'C', 'G']),
        (b'N', vec!['A', 'C', 'G', 'T']),
    ]
    .iter()
    .cloned()
    .collect();

    let mut all_kmers: HashMap<u128, HashMap<u128, u32>> = HashMap::new();
    let mut all_indexes: HashMap<String, u32> = HashMap::new();
    let mut index_map: HashMap<u32, String> = HashMap::new();

    let kmer_iter = ska_array.iter();

    for (int_kmer, int_middle_base_vec) in kmer_iter {
        // select non-ubiquitous split k-mers (ie, absent in at least one sample)
        if int_middle_base_vec.contains(&45) {
            let (kmer_left, kmer_right) = decode_kmer(len_kmer, int_kmer, mask, mask);

            // combine samples by middle-base by using degenerate code
            let mut middle_2_samples: HashMap<char, Vec<String>> = HashMap::new();
            for (i, nucl) in int_middle_base_vec.iter().enumerate() {
                // if middle-base not absent
                if *nucl != 45 {
                    for &new_nucl in &degenerate_code[nucl] {
                        middle_2_samples
                            .entry(new_nucl)
                            .or_default()
                            .push(i.to_string());
                    }
                }
            }

            for (nucl, l_indexes) in middle_2_samples.iter() {
                let str_sample_id = l_indexes.join("|");
                let value_index: u32;

                // get index to save species list (as String)
                if let Some(index) = all_indexes.get(&str_sample_id) {
                    // Key exists, use the existing index
                    value_index = *index;
                } else {
                    // Key does not exist, insert a new entry with the next index
                    value_index = all_indexes.len() as u32;
                    all_indexes.insert(str_sample_id.to_string(), value_index.clone());
                    index_map.insert(value_index.clone(), str_sample_id.to_string());
                }

                // save kmers in coloured de Bruijn graph
                let full_kmer = format!("{}{}{}", kmer_left, nucl, kmer_right);

                let encoded_kmer_1 = encode_kmer(&full_kmer[..full_kmer.len() - 1].to_string());
                let encoded_kmer_2 = encode_kmer(&full_kmer[1..].to_string());

                // uncomment to print network
                // println!("{}	{}	", &full_kmer[..full_kmer.len() - 1].to_string(), &full_kmer[1..].to_string());

                all_kmers
                    .entry(encoded_kmer_1)
                    .or_default()
                    .insert(encoded_kmer_2, value_index.clone());

                let rc_kmer = rev_compl(&full_kmer);
                let rc_encoded_kmer_1 = encode_kmer(&rc_kmer[..full_kmer.len() - 1].to_string());
                let rc_encoded_kmer_2 = encode_kmer(&rc_kmer[1..].to_string());

                // uncomment to print network
                // println!("{}	{}	", &rc_kmer[..full_kmer.len() - 1].to_string(), &rc_kmer[1..].to_string());

                all_kmers
                    .entry(rc_encoded_kmer_1)
                    .or_default()
                    .insert(rc_encoded_kmer_2, value_index);
            }
        }
    }

    // sanity check: test if the number of sample combinations is too high to be stored as u64 integers
    let max_u32_value: usize = u32::MAX as usize;
    if all_indexes.len() > max_u32_value {
        eprintln!(
            "\nError: the number of sample combinations is too high to be stored as u64 integers\n"
        );
        std::process::exit(1);
    }

    println!("     . {} nodes (includes rev-compl)", all_kmers.len());
    println!("     . {} sample combinations", all_indexes.len());
    (len_kmer, sample_names, all_kmers, index_map)
}
