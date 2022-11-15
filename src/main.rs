
use rayon::prelude::*;

pub mod ska_dict;
use crate::ska_dict::{SkaDict, MergeSkaDict};
use std::time::Instant;

fn main() {
    let n_threads = 4;
    rayon::ThreadPoolBuilder::new().num_threads(n_threads).build_global().unwrap();
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
    let kmer_size: u8 = 31;
    let rc = true;

    let start = Instant::now();
    let mut ska_dicts: Vec<SkaDict> = Vec::new();
    for file_it in file_list.iter().enumerate() {
        let (idx, (name, filename)) = file_it;
        ska_dicts.push(SkaDict::new(idx, filename, name, rc))
    }
    let build = Instant::now();

    let mut merged_dict = MergeSkaDict::new(ska_dicts.len());
    for ska_dict in &mut ska_dicts {
        merged_dict.append(ska_dict);
    }
    let merge = Instant::now();

    let ska_dict_p =
         file_list
         .par_iter()
         .enumerate()
         .map(|(idx, (name, filename))| SkaDict::new(idx, filename, name, rc))
         .fold(|| MergeSkaDict::new(file_list.len()),
                |mut a: MergeSkaDict, b: SkaDict| {a.append(&b); a})
         .reduce(|| MergeSkaDict::new(file_list.len()),
                  |mut a: MergeSkaDict, mut b: MergeSkaDict| { a.merge(&mut b); a });
    let parallel = Instant::now();

    print!("{}", merged_dict);
    print!("{}", ska_dict_p);
    println!("build:\t{}ms\nmerge:\t{}ms\npara:\t{}ms",
             build.duration_since(start).as_millis(),
             merge.duration_since(build).as_millis(),
             parallel.duration_since(merge).as_millis());
}
