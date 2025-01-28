use bit_set::BitSet;
use hashbrown::HashMap;

use dashmap::DashMap;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;

use crate::merge_ska_array::MergeSkaArray;
use crate::ska_dict::bit_encoding::{decode_kmer, UInt};

type KmerGraph<IntT> = HashMap<IntT, Vec<IntT>>;
type KmerSamples<IntT> = HashMap<IntT, BitSet>;

pub fn build_graph<IntT: for<'a> UInt<'a>>(
    ska_array: MergeSkaArray<IntT>,
    nb_threads: usize,
) -> (usize, Vec<String>, KmerGraph<IntT>, KmerSamples<IntT>) {
    // let arguments = CONFIG.get().unwrap();

    // log::info!(" # read file {}", config.input_file);

    // read the skf file and load split-kmers (ska_array), kmer length and sample names
    // let ska_array = load_array::<IntT>(&[arguments.input_file.to_string()], arguments.nb_threads)
    //     .expect("\nerror: could not read the skf file\n\n");
    let sample_names = ska_array.names().to_vec();
    let len_kmer = ska_array.kmer_len();
    let mask: IntT = IntT::skalo_mask(len_kmer);

    log::info!("     . {}-mers", len_kmer);
    log::info!("     . {} samples", sample_names.len());

    log::info!(" # build colored de Bruijn graph");

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

    let all_kmers: DashMap<IntT, Vec<IntT>> = DashMap::new();
    let kmer_samples: DashMap<IntT, BitSet> = DashMap::new();

    ThreadPoolBuilder::new()
        .num_threads(nb_threads)
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

                let encoded_kmer_1 = IntT::encode_kmer_str(&full_kmer[..len_kmer - 1]);
                let encoded_kmer_2 = IntT::encode_kmer_str(&full_kmer[1..]);

                all_kmers
                    .entry(encoded_kmer_1)
                    .or_default()
                    .push(encoded_kmer_2);

                all_kmers
                    .entry(IntT::rev_comp(encoded_kmer_2, len_kmer - 1))
                    .or_default()
                    .push(IntT::rev_comp(encoded_kmer_1, len_kmer - 1));

                let encode_full = IntT::encode_kmer_str(&full_kmer);
                kmer_samples
                    .entry(encode_full)
                    .or_insert_with(|| bitset_samples.clone());
                kmer_samples
                    .entry(IntT::rev_comp(encode_full, len_kmer))
                    .or_insert_with(|| bitset_samples.clone());
            }
        });

    let all_kmers: KmerGraph<IntT> = all_kmers.into_iter().collect();
    let kmer_samples: KmerSamples<IntT> = kmer_samples.into_iter().collect();

    log::info!("     . {} nodes", all_kmers.len());

    (len_kmer, sample_names, all_kmers, kmer_samples)
}
