
use rayon::prelude::*;

use std::fs::File;
use std::io::{BufWriter, Write};

pub mod ska_dict;
use crate::ska_dict::SkaDict;

pub mod merge_ska_dict;
use crate::merge_ska_dict::MergeSkaDict;

pub mod merge_ska_array;
use crate::merge_ska_array::MergeSkaArray;

use std::time::Instant;

fn main() {
    let n_threads = 4;
    rayon::ThreadPoolBuilder::new().num_threads(n_threads).build_global().unwrap();
    let small_file = vec![
        ("sample1", "test1.fa"),
        ("sample2", "test2.fa")
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

    let start = Instant::now();
    let mut ska_dicts: Vec<SkaDict> = Vec::new();
    for file_it in small_file.iter().enumerate() {
        let (idx, (name, filename)) = file_it;
        ska_dicts.push(SkaDict::new(kmer_size, idx, filename, name, rc))
    }
    let build = Instant::now();

    let mut merged_dict = MergeSkaDict::new(kmer_size, ska_dicts.len());
    for ska_dict in &mut ska_dicts {
        merged_dict.append(ska_dict);
    }
    let merge = Instant::now();

    let ska_dict_p =
    small_file
         .par_iter()
         .enumerate()
         .map(|(idx, (name, filename))| SkaDict::new(kmer_size, idx, filename, name, rc))
         .fold(|| MergeSkaDict::new(kmer_size, small_file.len()),
                |mut a: MergeSkaDict, b: SkaDict| {a.append(&b); a})
         .reduce(|| MergeSkaDict::new(kmer_size, small_file.len()),
                  |mut a: MergeSkaDict, mut b: MergeSkaDict| { a.merge(&mut b); a });
    let parallel = Instant::now();

    //print!("{}", merged_dict);
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
    ska_array.filter(ska_array.nsamples(), true);
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
}

// Command line
// Input commands
// ska build <seq.fa> --compress
// ska build -l rfile.txt --compress

// Manipulation commands
// ska merge file.skf...
// ska delete file.skf <names>...
// ska delete file.skf -l name_list.txt
// ska weed file.skf --weed-seqs list.fa

// Ouput commands (these take either one skf file, or an rfile, or a list on the command line)
// ska align <seq.fa>... --save-skf <-o output.aln> <--min 0.9>
// ska align -l rfile.txt --save-skf <-o output.aln> <--min 0.9>
// ska map ref.fa <seq.fa>... --save-skf <-o output.aln> <--min 0.9>
// ska map -l rfile.txt --save-skf <-o output.aln> <--min 0.9>

// Printing/debug commands
// ska nk file.skf // prints the number of k-mers (he he)
// ska info file.skf <--bases> // prints the k-mer size, number of samples, names, number of k-mers <number of patterns, prop A/C/G/T, GC-bias>
// ska humanise file.skf // prints the split k-mers after decoding them