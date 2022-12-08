use std::fs::File;
use std::io::{stdout, BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::time::Instant;

use simple_logger;

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
            .enumerate()
            .map(|(idx, (name, filename))| SkaDict::new(k, idx, filename, name, rc))
            .collect();
    } else {
        for file_it in input_files.iter().enumerate() {
            let (idx, (name, filename)) = file_it;
            ska_dicts.push(SkaDict::new(k, idx, filename, name, rc))
        }
    }

    // Merge indexes
    log::debug!("Merging skf dicts");
    let mut merged_dict = MergeSkaDict::new(k, ska_dicts.len(), rc);
    for ska_dict in &mut ska_dicts {
        merged_dict.append(ska_dict);
    }
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
        let merged_dict = build_and_merge(&input_files, DEFAULT_KMER, DEFAULT_STRAND, threads);
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
    simple_logger::init_with_env().unwrap();
    let args = cli_args();

    let start = Instant::now();
    match &args.command {
        Commands::Build {
            seq_files,
            file_list,
            output,
            k,
            single_strand,
            threads,
        } => {
            let rc = !*single_strand;
            // Read input
            let input_files = get_input_list(file_list, seq_files);

            // Build, merge
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
            log::debug!("Applying filters");
            let filter_threshold = f64::ceil(ska_array.nsamples() as f64 * *min_freq) as usize;
            ska_array.filter(filter_threshold, *const_sites);

            // Write out to file/stdout
            let mut out_stream = set_ostream(output);
            log::debug!("Writing alignment");
            write!(&mut out_stream, "{}", ska_array).unwrap();
        }
        Commands::Map {
            reference,
            input,
            output,
            output_format,
            threads,
        } => {
            log::debug!("Loading skf as dictionary");
            let ska_dict = load_array(input, *threads).to_dict();

            log::debug!("Making skf of reference");
            let mut ska_ref = RefSka::new(ska_dict.ksize(), &reference, ska_dict.rc());

            log::debug!("Mapping");
            ska_ref.map(&ska_dict);

            let mut out_stream = set_ostream(output);
            match output_format {
                FileType::Aln => {
                    log::debug!("Writing alignment");
                    write!(&mut out_stream, "{}", ska_ref).unwrap();
                }
                FileType::Bcf => {
                    log::debug!("Writing BCF");
                    ska_ref.write_vcf(&mut out_stream, true);
                }
                FileType::Vcf => {
                    log::debug!("Writing VCF");
                    ska_ref.write_vcf(&mut out_stream, false);
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
        } => {}
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
    /*
    let small_file = vec![
        ("sample1", "N_test_1.fa"),
        ("sample2", "N_test_2.fa")
    ];
    let file_list = vec![
        ("BR1076_4336457", "assemblies/BR1076_4336457.contigs.fa"),
        ("NP7078_4383169", "assemblies/NP7078_4383169.contigs.fa"),
        ("LE4079_4336554", "assemblies/LE4079_4336554.contigs.fa"),
        ("PT8092_4383213", "assemblies/PT8092_4383213.contigs.fa"),
        ("ND6110_4382925", "assemblies/ND6110_4382925.contigs.fa"),
        ("3053_3610881", "assemblies/3053_3610881.contigs.fa"),
        ("3060_3610889", "assemblies/3060_3610889.contigs.fa"),
        ("GL3032_4336520", "assemblies/GL3032_4336520.contigs.fa"),
        ("PT8081_4383208", "assemblies/PT8081_4383208.contigs.fa"),
        ("NP7028_4382961", "assemblies/NP7028_4382961.contigs.fa"),
        ("093209_3736979", "assemblies/093209_3736979.contigs.fa")
    ];
    let kmer_size: usize = 31;
    let rc = true;
    let const_sites = false;

    let start = Instant::now();
    let mut ska_dicts: Vec<SkaDict> = Vec::new();
    for file_it in small_file.iter().enumerate() {
        let (idx, (name, filename)) = file_it;
        ska_dicts.push(SkaDict::new(kmer_size, idx, filename, name, rc))
    }
    let build = Instant::now();

    let mut merged_dict = MergeSkaDict::new(kmer_size, ska_dicts.len(), rc);
    for ska_dict in &mut ska_dicts {
        merged_dict.append(ska_dict);
    }
    let merge = Instant::now();

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

    print!("{}", merged_dict);
    // print!("{}", ska_dict_p);
    println!("build:\t{}ms\nmerge:\t{}ms\npara:\t{}ms",
             build.duration_since(start).as_millis(),
             merge.duration_since(build).as_millis(),
             parallel.duration_since(merge).as_millis());

    let convert = Instant::now();
    let mut ska_array = MergeSkaArray::new(&merged_dict);
    let deconvert = Instant::now();
    let _ska_dict_c = ska_array.to_dict();
    let deconvert_end = Instant::now();

    println!("encode:\t{}ms\ndecode:\t{}ms",
    deconvert.duration_since(convert).as_millis(),
    deconvert_end.duration_since(deconvert).as_millis());

    let io_start = Instant::now();
    let mut file = BufWriter::new(File::create("test.aln").unwrap());
    ska_array.filter(ska_array.nsamples(), const_sites);
    write!(&mut file, "{}", ska_array).unwrap();

    let save = Instant::now();
    ska_array.save("test.skf").unwrap();

    let load = Instant::now();
    let _ska_array_load = MergeSkaArray::load("test.skf").unwrap();
    let load_end = Instant::now();

    println!("write:\t{}ms\nsave:\t{}ms\nload:\t{}ms",
        save.duration_since(io_start).as_millis(),
        load.duration_since(save).as_millis(),
        load_end.duration_since(load).as_millis());
    */
}
