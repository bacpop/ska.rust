use hashbrown::{HashMap, HashSet};
//use std::time::Instant;
use bit_set::BitSet;

use crate::skalo::utils::{rev_compl_u128, DATA_INFO};

pub fn identify_good_kmers(
    all_kmers: &HashMap<u128, Vec<u128>>,
    kmer_2_samples: &HashMap<u128, BitSet>,
) -> (HashSet<u128>, HashSet<u128>) {
    println!(" # identify bubble extremities");

    let data_info = DATA_INFO.get().unwrap();

    let mut start_kmers: HashSet<u128> = HashSet::new();
    let mut end_kmers: HashSet<u128> = HashSet::new();

    //let start = Instant::now();

    // iterate over all_kmers
    for (kmer, next_kmers) in all_kmers.iter() {
        if next_kmers.len() > 1 {
            'i_loop: for (i, &kmer1) in next_kmers.iter().enumerate() {
                for &kmer2 in next_kmers.iter().skip(i + 1) {
                    let full_kmer1 = combine_kmers(*kmer, kmer1);
                    let full_kmer2 = combine_kmers(*kmer, kmer2);

                    let samples1 = kmer_2_samples.get(&full_kmer1).unwrap();
                    let samples2 = kmer_2_samples.get(&full_kmer2).unwrap();

                    if compare_samples(samples1, samples2) {
                        start_kmers.insert(*kmer);
                        end_kmers.insert(rev_compl_u128(*kmer, data_info.k_graph));

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
        eprintln!("\n      Error: there is no entry node in this graph, hence no variant.\n");
        std::process::exit(1);
    }

    println!("     . {} entry nodes", start_kmers.len());

    //let duration = start.elapsed();
    //println!("time taken: {:?}", duration);

    (start_kmers, end_kmers)
}

pub fn combine_kmers(encoded_kmer1: u128, encoded_kmer2: u128) -> u128 {
    // define the bit mask for extracting the last nucleotide of the second k-mer
    let last_nucleotide_mask: u128 = 0b11; // Mask for 2 bits

    // shift the first k-mer left by 2 bits to make space for the new nucleotide
    let shifted_kmer1 = encoded_kmer1 << 2;

    // extract the last nucleotide from the second k-mer
    let last_nucleotide = encoded_kmer2 & last_nucleotide_mask;

    // combine the two k-mers into a (k+1)-mer encoding
    shifted_kmer1 | last_nucleotide
}

fn compare_samples(bitset_1: &BitSet, bitset_2: &BitSet) -> bool {
    let diff_count = bitset_1.symmetric_difference(bitset_2).count();
    diff_count != 0
}
