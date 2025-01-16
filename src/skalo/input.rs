use hashbrown::HashMap;
//use std::time::Instant;
use bit_set::BitSet;

use dashmap::DashMap;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;

use crate::io_utils::load_array;
use crate::ska_dict::bit_encoding::decode_kmer;

use crate::skalo::utils::{encode_kmer, rev_compl_u128, CONFIG};

type KmerGraph = HashMap<u128, Vec<u128>>;
type KmerSamples = HashMap<u128, BitSet>;

pub fn read_input_file() -> (usize, Vec<String>, KmerGraph, KmerSamples) {
    let arguments = CONFIG.get().unwrap();

    println!(" # read file {}", arguments.input_file);

    //let start = Instant::now();

    // read the skf file and load split-kmers (ska_array), kmer length and sample names
    let ska_array = load_array::<u128>(&[arguments.input_file.to_string()], arguments.nb_threads)
        .expect("\nerror: could not read the skf file\n\n");
    let sample_names = ska_array.names().to_vec();
    let len_kmer = ska_array.kmer_len();
    let mask = (1 << (len_kmer * 2)) - 1;

    println!("     . {}-mers", len_kmer);
    println!("     . {} samples", sample_names.len());

    //let duration = start.elapsed();
    //println!("time taken: {:?}", duration);

    println!(" # build colored de Bruijn graph");

    //let start = Instant::now();

    // build De Bruijn graph
    let degenerate_code: HashMap<u8, Vec<char>> = [
        (b'A', vec!['A']),
        (b'T', vec!['T']),
        (b'G', vec!['G']),
        (b'C', vec!['C']),
        (b'M', vec!['A', 'C']),
        (b'S', vec!['C', 'G']),
        (b'W', vec!['A', 'T']),
        (b'R', vec!['A', 'G']),
        (b'Y', vec!['C', 'T']),
        (b'K', vec!['G', 'T']),
        (b'B', vec!['C', 'G', 'T']),
        (b'D', vec!['A', 'G', 'T']),
        (b'H', vec!['A', 'C', 'T']),
        (b'V', vec!['A', 'C', 'G']),
        (b'N', vec!['A', 'C', 'G', 'T']),
    ]
    .iter()
    .cloned()
    .collect();

    let all_kmers: DashMap<u128, Vec<u128>> = DashMap::new();
    let kmer_samples: DashMap<u128, BitSet> = DashMap::new();

    ThreadPoolBuilder::new()
        .num_threads(arguments.nb_threads)
        .build_global()
        .expect("failed to build the thread pool");

    let kmer_iter = ska_array.iter();

    kmer_iter
        .par_bridge()
        .for_each(|(int_kmer, int_middle_base_vec)| {
            let (kmer_left, kmer_right) = decode_kmer(len_kmer, int_kmer, mask, mask);

            // combine samples by middle-base using degenerate code
            let mut middle_2_samples: HashMap<char, Vec<u16>> = HashMap::with_capacity(4);
            for (i, nucl) in int_middle_base_vec.iter().enumerate() {
                if *nucl != 45 {
                    for &new_nucl in &degenerate_code[nucl] {
                        middle_2_samples
                            .entry(new_nucl)
                            .or_insert_with(|| Vec::with_capacity(sample_names.len()))
                            .push(i as u16);
                    }
                }
            }

            // build k-mers and save them to the graph
            let mut bitset_samples = BitSet::with_capacity(sample_names.len());
            let mut full_kmer = String::with_capacity(len_kmer);

            for (nucl, vec_indexes) in middle_2_samples.iter() {
                // crate BitSet of samples
                bitset_samples.clear();
                bitset_samples.extend(vec_indexes.iter().map(|&x| x as usize));

                // save k-mers (k-1 in graph and k for samples)
                full_kmer.clear();
                full_kmer.push_str(&kmer_left);
                full_kmer.push(*nucl);
                full_kmer.push_str(&kmer_right);

                let encoded_kmer_1 = encode_kmer(&full_kmer[..len_kmer - 1]);
                let encoded_kmer_2 = encode_kmer(&full_kmer[1..]);

                all_kmers
                    .entry(encoded_kmer_1)
                    .or_default()
                    .push(encoded_kmer_2);

                all_kmers
                    .entry(rev_compl_u128(encoded_kmer_2, len_kmer - 1))
                    .or_default()
                    .push(rev_compl_u128(encoded_kmer_1, len_kmer - 1));

                let encode_full = encode_kmer(&full_kmer);
                kmer_samples
                    .entry(encode_full)
                    .or_insert_with(|| bitset_samples.clone());
                kmer_samples
                    .entry(rev_compl_u128(encode_full, len_kmer))
                    .or_insert_with(|| bitset_samples.clone());
            }
        });

    let all_kmers: KmerGraph = all_kmers.into_iter().collect();
    let kmer_samples: KmerSamples = kmer_samples.into_iter().collect();

    //let duration = start.elapsed();
    //println!("time taken: {:?}", duration);
    println!("     . {} nodes", all_kmers.len());

    (len_kmer, sample_names, all_kmers, kmer_samples)
}
