use std::time::Instant;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

use regex::Regex;
use rayon::prelude::*;

pub mod ska_dict;
use crate::ska_dict::SkaDict;

pub mod merge_ska_dict;
use crate::merge_ska_dict::MergeSkaDict;

pub mod ska_ref;
pub mod merge_ska_array;
use crate::merge_ska_array::MergeSkaArray;

pub mod cli;
use crate::cli::{cli_args, Commands};

fn main() {
    let args = cli_args();

    let start = Instant::now();
    match &args.command {
        Commands::Build { seq_files, file_list, output, k, single_strand, threads } => {
            // Read input
            let rc = !single_strand;
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
                    let re = Regex::new(r"^(.+)\.(?i:fa|fasta)$").unwrap();
                    for file in seq_files.as_ref().unwrap() {
                        let caps = re.captures(file);
                        let name = match caps {
                            Some(capture) => capture[1].to_string(),
                            None => file.to_string()
                        };
                        input_files.push((name, file.to_string()));
                    }
                }
            }

            // Build indexes
            let mut ska_dicts: Vec<SkaDict> = Vec::new();
            ska_dicts.reserve(input_files.len());
            if *threads > 1 {
                rayon::ThreadPoolBuilder::new().num_threads(*threads).build_global().unwrap();
                ska_dicts =
                    input_files
                     .par_iter()
                     .enumerate()
                     .map(|(idx, (name, filename))| SkaDict::new(*k, idx, filename, name, rc))
                     .collect();
            } else {
                for file_it in input_files.iter().enumerate() {
                    let (idx, (name, filename)) = file_it;
                    ska_dicts.push(SkaDict::new(*k, idx, filename, name, rc))
                }
            }

            // Merge indexes
            let mut merged_dict = MergeSkaDict::new(*k, ska_dicts.len(), rc);
            for ska_dict in &mut ska_dicts {
                merged_dict.append(ska_dict);
            }

            // Save
            let ska_array = MergeSkaArray::new(&merged_dict);
            ska_array.save(format!("{output}.skf").as_str()).expect("Failed to save output file");
        }
        Commands::Nk { skf_file, full_info } => {
            let ska_array_load = MergeSkaArray::load(skf_file).unwrap();
            let ska_dict = ska_array_load.to_dict();
            println!("{}", ska_dict);
            if *full_info {
                println!("{:?}", ska_dict);
            }
        }
        _ => {}
    }
    let end = Instant::now();

    println!("SKA done in {}ms",
        end.duration_since(start).as_millis());
    println!("⬛⬜⬛⬜⬛⬜⬛⬜⬛⬜");
    println!("⬜⬛⬜⬛⬜⬛⬜⬛⬜⬛");
    println!("⬛⬜⬛⬜⬛⬜⬛⬜⬛⬜");
    println!("⬜⬛⬜⬛⬜⬛⬜⬛⬜⬛");
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

