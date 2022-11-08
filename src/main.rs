
use rayon::prelude::*;

pub mod ska_dict;
use crate::ska_dict::SkaDict;
use std::time::Instant;

// A better way to do this would be to start with some classes
// ska_dict: HashMap<u64, u8>, with read method, to multi_ska
// multi_ska: HashMap<u64, <bitset, bitset, bitset, bitset>>
// methods: merge, merge (w/ ska_dict), filter, write

// Then we can have
// use rayon::prelude::*;

// let ska_dict =
//     file_list
//     .par_iter()
//     .map(|f: &str| SkaDict.from(f))
//     .reduce(|| SkaDict::new(),
//             |mut a: SkaDict, b: SkaDict| { a.merge(&b); a });

// Easiest to start with bitvec of length 1
// Then on merge
// if match: concat the two bitvecs
// if no match: concat found with zeros * n_samples in no match

fn main() {
    let start = Instant::now();
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

    //let ska_file = SkaDict::new(filename, "19183_4#48", true);
    //print!("{}", ska_file);
    let ska_dict =
        file_list
        .par_iter()
        .map(|(name, filename)| SkaDict::new(filename, name, rc))
        .reduce(|| SkaDict::default(),
                 |mut a: SkaDict, mut b: SkaDict| { a.merge(&mut b); a });
    print!("{}", ska_dict);
    let end = Instant::now();
    println!("time taken: {}ms",
             end.duration_since(start).as_millis());
}
