use std::fs::File;
use std::io::{stdout, BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::time::Instant;

use simple_logger;
use indicatif::{ProgressBar, ProgressIterator, ParallelProgressIterator};

use rayon::prelude::*;
use regex::Regex;

pub mod ska_dict;
use crate::ska_dict::SkaDict;

pub mod merge_ska_dict;
use crate::merge_ska_dict::MergeSkaDict;

pub mod ska_ref;
use crate::ska_ref::RefSka;
pub mod merge_ska_array;
use crate::merge_ska_array::MergeSkaArray;

pub mod cli;
use crate::cli::*;

fn read_input_fastas(seq_files: &Vec<String>) -> Vec<(String, String)> {
    let mut input_files = Vec::new();
    let re = Regex::new(r"^(.+)\.(?i:fa|fasta|fastq|fastq\.gz)$").unwrap();
    for file in seq_files {
        let caps = re.captures(file);
        let name = match caps {
            Some(capture) => capture[1].to_string(),
            None => file.to_string(),
        };
        input_files.push((name, file.to_string()));
    }
    return input_files;
}

fn build_and_merge(
    input_files: &Vec<(String, String)>,
    k: usize,
    rc: bool,
    threads: usize,
) -> MergeSkaDict {
    // Build indexes
    log::debug!("Building skf dicts from sequence input");
    let mut ska_dicts: Vec<SkaDict> = Vec::new();
    ska_dicts.reserve(input_files.len());
    if threads > 1 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .unwrap();
        ska_dicts = input_files
            .par_iter()
            .progress_count(input_files.len() as u64)
            .enumerate()
            .map(|(idx, (name, filename))| SkaDict::new(k, idx, filename, name, rc))
            .collect();
    } else {
        for file_it in input_files.iter().progress().enumerate() {
            let (idx, (name, filename)) = file_it;
            ska_dicts.push(SkaDict::new(k, idx, filename, name, rc))
        }
    }

    // Parallel build & merge is slower
    /*
    let n_threads = 4;
    rayon::ThreadPoolBuilder::new().num_threads(n_threads).build_global().unwrap();
    let ska_dict_p =
    small_file
         .par_iter()
         .enumerate()
         .map(|(idx, (name, filename))| SkaDict::new(kmer_size, idx, filename, name, rc))
         .fold(|| MergeSkaDict::new(kmer_size, small_file.len(), rc),
                |mut a: MergeSkaDict, b: SkaDict| {a.append(&b); a})
         .reduce(|| MergeSkaDict::new(kmer_size, small_file.len(), rc),
                  |mut a: MergeSkaDict, mut b: MergeSkaDict| { a.merge(&mut b); a });
    let parallel = Instant::now();
    */

    // Merge indexes
    log::debug!("Merging skf dicts");
    let mut merged_dict = MergeSkaDict::new(k, ska_dicts.len(), rc);
    let bar = ProgressBar::new(ska_dicts.len() as u64);
    for ska_dict in &mut ska_dicts {
        merged_dict.append(ska_dict);
        bar.inc(1);
    }
    bar.finish();
    return merged_dict;
}

fn load_array(input: &Vec<String>, threads: usize) -> MergeSkaArray {
    // Obtain a merged ska array
    let ska_array: MergeSkaArray;
    if input.len() == 1 {
        log::debug!("Single file as input, trying to load as skf");
        ska_array = MergeSkaArray::load(input[0].as_str()).unwrap();
    } else {
        log::debug!("Multiple files as input, running ska build with default settings");
        let input_files = read_input_fastas(input);
        let merged_dict = build_and_merge(&input_files, DEFAULT_KMER, !DEFAULT_STRAND, threads);
        ska_array = MergeSkaArray::new(&merged_dict);
    }
    return ska_array;
}

// Write out to file/stdout
fn set_ostream(oprefix: &Option<String>) -> BufWriter<Box<dyn Write>> {
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

fn get_input_list(
    file_list: &Option<String>,
    seq_files: &Option<Vec<String>>,
) -> Vec<(String, String)> {
    // Read input
    let mut input_files: Vec<(String, String)> = Vec::new();
    match file_list {
        Some(files) => {
            let f = File::open(files).expect("Unable to open file_list");
            let f = BufReader::new(f);
            for line in f.lines() {
                let line = line.expect("Unable to read line in file_list");
                let fields: Vec<&str> = line.split_whitespace().collect();
                input_files.push((fields[0].to_string(), fields[1].to_string()));
            }
        }
        None => {
            input_files = read_input_fastas(seq_files.as_ref().unwrap());
        }
    }
    return input_files;
}

fn main() {
    log::debug!("Loading skf as dictionary");
    let ska_dict = load_array(&vec!["test_1.fa".to_string(), "test_2.fa".to_string()], 1).to_dict();

    log::debug!("Making skf of reference k={} rc={}", ska_dict.kmer_len(), ska_dict.rc());
    let mut ska_ref = RefSka::new(ska_dict.kmer_len(), "test_ref.fa", ska_dict.rc());

    log::debug!("Mapping");
    ska_ref.map(&ska_dict);

    let mut out_stream = set_ostream(&None);
    ska_ref.write_aln(&mut out_stream).expect("Failed to write output alignment");
    ska_ref.write_vcf(&mut out_stream).expect("Failed to write output VCF");
}

fn main2() {
    let args = cli_args();
    if args.verbose {
        simple_logger::init_with_level(log::Level::Debug).unwrap();
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
            threads
        } => {
            // Read input
            let input_files = get_input_list(file_list, seq_files);

            // Build, merge
            let rc = !*single_strand;
            let merged_dict = build_and_merge(&input_files, *k, rc, *threads);

            // Save
            log::debug!("Converting to array representation and saving");
            let ska_array = MergeSkaArray::new(&merged_dict);
            ska_array
                .save(format!("{output}.skf").as_str())
                .expect("Failed to save output file");
        }
        Commands::Align {
            input,
            output,
            min_freq,
            const_sites,
            threads,
        } => {
            let mut ska_array = load_array(input, *threads);

            // Apply filters
            let filter_threshold = f64::ceil(ska_array.nsamples() as f64 * *min_freq) as usize;
            log::debug!("Applying filters: threshold={filter_threshold} const_sites={const_sites}");
            ska_array.filter(filter_threshold, *const_sites);

            // Write out to file/stdout
            let mut out_stream = set_ostream(output);
            log::debug!("Writing alignment");
            ska_array.write_fasta(&mut out_stream).expect("Couldn't write output fasta");
        }
        Commands::Map {
            reference,
            input,
            output,
            format,
            threads,
        } => {
            log::debug!("Loading skf as dictionary");
            let ska_dict = load_array(input, *threads).to_dict();

            log::debug!("Making skf of reference k={} rc={}", ska_dict.kmer_len(), ska_dict.rc());
            let mut ska_ref = RefSka::new(ska_dict.kmer_len(), &reference, ska_dict.rc());

            log::debug!("Mapping");
            ska_ref.map(&ska_dict);

            let mut out_stream = set_ostream(output);
            match format {
                FileType::Aln => {
                    log::debug!("Writing alignment");
                    ska_ref.write_aln(&mut out_stream).expect("Failed to write output alignment");
                }
                FileType::Vcf => {
                    log::debug!("Writing VCF");
                    ska_ref.write_vcf(&mut out_stream).expect("Failed to write output VCF");
                }
            }
        }
        Commands::Merge { skf_files, output } => {
            if skf_files.len() < 2 {
                panic!("Need at least two files to merge");
            }

            let first_array =
                MergeSkaArray::load(&skf_files[0]).expect("Failed to load input file");
            let mut merged_dict = first_array.to_dict();
            for file_idx in 1..skf_files.len() {
                log::debug!("Merging alignment {file_idx}");
                let next_array =
                    MergeSkaArray::load(&skf_files[file_idx]).expect("Failed to load input file");
                merged_dict.merge(&mut next_array.to_dict());
            }
            let merged_array = MergeSkaArray::new(&merged_dict);
            merged_array
                .save(format!("{output}.skf").as_str())
                .expect("Failed to save output file");
        }
        Commands::Delete {
            skf_file,
            file_list,
            names,
        } => {
            log::debug!("Loading skf file");
            let mut ska_array =
                MergeSkaArray::load(skf_file).expect("Could not load input skf file");
            let input_files = get_input_list(file_list, names);

            log::debug!("Deleting samples");
            ska_array.delete_samples(&input_files);

            log::debug!("Saving modified skf file");
            ska_array
                .save(skf_file)
                .expect("Could not save modified array");
        }
        Commands::Weed {
            skf_file,
            weed_file,
        } => {
            log::debug!("Loading skf as dictionary");
            let mut ska_dict = MergeSkaArray::load(skf_file.as_str()).unwrap().to_dict();

            log::debug!("Making skf of weed file k={} rc={}", ska_dict.kmer_len(), ska_dict.rc());
            let ska_weed = RefSka::new(ska_dict.kmer_len(), &weed_file, ska_dict.rc());

            log::debug!("Removing weed k-mers");
            ska_dict.weed(&ska_weed);

            log::debug!("Saving modified skf file");
            let ska_array = MergeSkaArray::new(&ska_dict);
            ska_array
                .save(skf_file.as_str())
                .expect("Failed to save output file");
        }
        Commands::Nk {
            skf_file,
            full_info,
        } => {
            log::debug!("Printing basic info");
            let ska_array_load = MergeSkaArray::load(skf_file).unwrap();
            let ska_dict = ska_array_load.to_dict();
            println!("{}", ska_dict);

            if *full_info {
                log::debug!("Printing full info");
                println!("{:?}", ska_dict);
            }
        }
    }
    let end = Instant::now();

    eprintln!("SKA done in {}ms", end.duration_since(start).as_millis());
    eprintln!("⬛⬜⬛⬜⬛⬜⬛⬜⬛⬜");
    eprintln!("⬜⬛⬜⬛⬜⬛⬜⬛⬜⬛");
    eprintln!("⬛⬜⬛⬜⬛⬜⬛⬜⬛⬜");
    eprintln!("⬜⬛⬜⬛⬜⬛⬜⬛⬜⬛");
    log::debug!("Complete");
}
