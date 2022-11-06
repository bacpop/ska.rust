
use rayon::prelude::*;

pub mod ska_dict;

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
    let filename = "19183_4#48.contigs_velvet.fa";
    let kmer_size: u8 = 31;

    let ska_file = SkaDict::new(filename, "19183_4#48", true);
    print!("{}", ska_file);
    // eventually
    // let ska_dict =
    //     file_list
    //     .par_iter()
    //     .map(|f: &str| SkaDict.from(f))
    //     .reduce(|| SkaDict::Default(),
    //             |mut a: SkaDict, b: SkaDict| { a.merge(&b); a });
}
