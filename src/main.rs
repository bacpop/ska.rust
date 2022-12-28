use std::time::Instant;

use simple_logger;

pub mod ska_dict;
pub mod merge_ska_dict;
use crate::merge_ska_dict::build_and_merge;

pub mod ska_ref;
use crate::ska_ref::RefSka;
pub mod merge_ska_array;
use crate::merge_ska_array::MergeSkaArray;

pub mod cli;
use crate::cli::*;

pub mod io_utils;
use crate::io_utils::*;

fn main2() {
    log::info!("Loading skf as dictionary");
    let ska_dict = load_array(&vec!["test_1.fa".to_string(), "test_2.fa".to_string()], 1).to_dict();

    log::info!("Making skf of reference k={} rc={}", ska_dict.kmer_len(), ska_dict.rc());
    let mut ska_ref = RefSka::new(ska_dict.kmer_len(), "test_ref.fa", ska_dict.rc());

    log::info!("Mapping");
    ska_ref.map(&ska_dict);

    let mut out_stream = set_ostream(&None);
    ska_ref.write_aln(&mut out_stream).expect("Failed to write output alignment");
    ska_ref.write_vcf(&mut out_stream).expect("Failed to write output VCF");
}

fn main() {
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
            threads
        } => {
            // Read input
            let input_files = get_input_list(file_list, seq_files);

            // Build, merge
            let rc = !*single_strand;
            let merged_dict = build_and_merge(&input_files, *k, rc, *min_count, *min_qual, *threads);


            // Save
            log::info!("Converting to array representation and saving");
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
            // In debug mode (cannot be set from CLI, give details)
            log::debug!("{}", format!("{}", ska_array));

            // Apply filters
            let filter_threshold = f64::ceil(ska_array.nsamples() as f64 * *min_freq) as usize;
            log::info!("Applying filters: threshold={filter_threshold} const_sites={const_sites}");
            ska_array.filter(filter_threshold, *const_sites);

            // Write out to file/stdout
            let mut out_stream = set_ostream(output);
            log::info!("Writing alignment");
            ska_array.write_fasta(&mut out_stream).expect("Couldn't write output fasta");
        }
        Commands::Map {
            reference,
            input,
            output,
            format,
            threads,
        } => {
            log::info!("Loading skf as dictionary");
            let ska_dict = load_array(input, *threads).to_dict();

            log::info!("Making skf of reference k={} rc={}", ska_dict.kmer_len(), ska_dict.rc());
            let mut ska_ref = RefSka::new(ska_dict.kmer_len(), &reference, ska_dict.rc());

            log::info!("Mapping");
            ska_ref.map(&ska_dict);

            let mut out_stream = set_ostream(output);
            match format {
                FileType::Aln => {
                    log::info!("Writing alignment");
                    ska_ref.write_aln(&mut out_stream).expect("Failed to write output alignment");
                }
                FileType::Vcf => {
                    log::info!("Writing VCF");
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
                log::info!("Merging alignment {file_idx}");
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
            log::info!("Loading skf file");
            let mut ska_array =
                MergeSkaArray::load(skf_file).expect("Could not load input skf file");
            let input_files = get_input_list(file_list, names);

            log::info!("Deleting samples");
            let input_names = input_files.iter().map(|t| t.0.to_owned()).collect();
            ska_array.delete_samples(&input_names);

            log::info!("Saving modified skf file");
            ska_array
                .save(skf_file)
                .expect("Could not save modified array");
        }
        Commands::Weed {
            skf_file,
            weed_file,
        } => {
            log::info!("Loading skf file");
            let mut ska_array = MergeSkaArray::load(skf_file.as_str()).unwrap();

            log::info!("Making skf of weed file k={} rc={}", ska_array.kmer_len(), ska_array.rc());
            let ska_weed = RefSka::new(ska_array.kmer_len(), &weed_file, ska_array.rc());

            log::info!("Removing weed k-mers");
            ska_array.weed(&ska_weed);

            log::info!("Saving modified skf file");
            ska_array
                .save(skf_file.as_str())
                .expect("Failed to save output file");
        }
        Commands::Nk {
            skf_file,
            full_info,
        } => {
            log::info!("Printing basic info");
            let ska_array = MergeSkaArray::load(skf_file).unwrap();
            println!("{}", ska_array);

            if *full_info {
                log::info!("Printing full info");
                println!("{:?}", ska_array);
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
