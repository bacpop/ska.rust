use std::fs::File;
use std::io::{stdout, BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use regex::Regex;

use crate::merge_ska_array::MergeSkaArray;
use crate::merge_ska_dict::{build_and_merge, InputFastx};

use crate::cli::{DEFAULT_KMER, DEFAULT_MINCOUNT, DEFAULT_MINQUAL, DEFAULT_STRAND};

pub fn read_input_fastas(seq_files: &[String]) -> Vec<InputFastx> {
    let mut input_files = Vec::new();
    let re = Regex::new(r"^(.+)\.(?i:fa|fasta|fastq|fastq\.gz)$").unwrap();
    for file in seq_files {
        let caps = re.captures(file);
        let name = match caps {
            Some(capture) => capture[1].to_string(),
            None => file.to_string(),
        };
        input_files.push((name, file.to_string(), None));
    }
    return input_files;
}

pub fn load_array(input: &[String], threads: usize) -> MergeSkaArray {
    // Obtain a merged ska array
    let ska_array: MergeSkaArray;
    if input.len() == 1 {
        log::info!("Single file as input, trying to load as skf");
        ska_array = MergeSkaArray::load(input[0].as_str()).unwrap();
    } else {
        log::info!("Multiple files as input, running ska build with default settings");
        let input_files = read_input_fastas(input);
        let merged_dict = build_and_merge(
            &input_files,
            DEFAULT_KMER,
            !DEFAULT_STRAND,
            DEFAULT_MINCOUNT,
            DEFAULT_MINQUAL,
            threads,
        );
        ska_array = MergeSkaArray::new(&merged_dict);
    }
    return ska_array;
}

// Write out to file/stdout
pub fn set_ostream(oprefix: &Option<String>) -> BufWriter<Box<dyn Write>> {
    let out_writer = match oprefix {
        Some(prefix) => {
            let path = Path::new(prefix);
            Box::new(File::create(&path).unwrap()) as Box<dyn Write>
        }
        None => Box::new(stdout()) as Box<dyn Write>,
    };
    let out_stream = BufWriter::new(out_writer);
    return out_stream;
}

pub fn get_input_list(
    file_list: &Option<String>,
    seq_files: &Option<Vec<String>>,
) -> Vec<InputFastx> {
    // Read input
    let mut input_files: Vec<InputFastx> = Vec::new();
    match file_list {
        Some(files) => {
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
        }
        None => {
            input_files = read_input_fastas(seq_files.as_ref().unwrap());
        }
    }
    return input_files;
}
