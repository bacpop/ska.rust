use bit_set::BitSet;
use hashbrown::{HashMap, HashSet};

use crate::ska_dict::bit_encoding::UInt;
use crate::skalo::utils::DataInfo;

pub fn identify_good_kmers<IntT: for<'a> UInt<'a>>(
    all_kmers: &HashMap<IntT, Vec<IntT>>,
    kmer_2_samples: &HashMap<IntT, BitSet>,
    data_info: &DataInfo,
) -> (HashSet<IntT>, HashSet<IntT>) {
    log::info!(" # identify bubble extremities");

    let mut start_kmers: HashSet<IntT> = HashSet::new();
    let mut end_kmers: HashSet<IntT> = HashSet::new();

    // iterate over all_kmers
    for (kmer, next_kmers) in all_kmers.iter() {
        if next_kmers.len() > 1 {
            'i_loop: for (i, &kmer1) in next_kmers.iter().enumerate() {
                for &kmer2 in next_kmers.iter().skip(i + 1) {
                    let full_kmer1 = IntT::combine_kmers(*kmer, kmer1);
                    let full_kmer2 = IntT::combine_kmers(*kmer, kmer2);

                    let samples1 = kmer_2_samples.get(&full_kmer1).unwrap();
                    let samples2 = kmer_2_samples.get(&full_kmer2).unwrap();

                    if compare_samples(samples1, samples2) {
                        start_kmers.insert(*kmer);
                        end_kmers.insert(IntT::rev_comp(*kmer, data_info.k_graph));

                        //uncomment to print network
                        /*
                        use crate::utils::{rev_compl, decode_kmer};
                        let dna = decode_kmer(kmer.clone(), data_info.k_graph);
                        let rc = rev_compl(&dna);

                        println!("{}	{}	red	", &kmer, &kmer);
                        println!("{}	{}	red	", &rev_compl_u128(*kmer, data_info.k_graph), &rev_compl_u128(*kmer, data_info.k_graph));
                        */

                        break 'i_loop;
                    }
                }
            }
        }
    }

    // exit program if no extremity found (eg, cases of weeded skf files)
    if start_kmers.is_empty() {
        log::error!("\n      Error: there is no entry node in this graph, hence no variant.\n");
        std::process::exit(1);
    }

    log::info!("{} entry nodes", start_kmers.len());

    (start_kmers, end_kmers)
}

fn compare_samples(bitset_1: &BitSet, bitset_2: &BitSet) -> bool {
    let diff_count = bitset_1.symmetric_difference(bitset_2).count();
    diff_count != 0
}
