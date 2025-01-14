use crate::io_utils::load_array;
use crate::ska_dict::bit_encoding::decode_kmer;
use hashbrown::HashMap;

use crate::utils::{encode_kmer, rev_compl};

pub fn read_input_file(
    input_file: &str,
) -> (
    usize,
    Vec<String>,
    HashMap<u128, HashMap<u128, u32>>,
    HashMap<u32, String>,
) {
    println!(" # read file '{}'", input_file);

    // read the skf file and load split-kmers (ska_array), kmer length and sample names
    let ska_array = load_array::<u128>(&[input_file.to_string()], 1)
        .expect("\nerror: coud not read the skf file\n\n");
    let sample_names = ska_array.names().to_vec();
    let len_kmer = ska_array.kmer_len();
    let mask = (1 << (len_kmer * 2)) - 1;

    println!("     . {}-mers", len_kmer);
    println!("     . {} samples", sample_names.len());

    println!(" # build colored de Bruijn graph");

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

    let mut all_kmers: HashMap<u128, HashMap<u128, u32>> = HashMap::new();
    let mut all_indexes: HashMap<String, u32> = HashMap::new();
    let mut index_map: HashMap<u32, String> = HashMap::new();

    let kmer_iter = ska_array.iter();

    for (int_kmer, int_middle_base_vec) in kmer_iter {
        // select non-ubiquitous split k-mers (ie, absent in at least one sample)
        if int_middle_base_vec.contains(&45) {
            let (kmer_left, kmer_right) = decode_kmer(len_kmer, int_kmer, mask, mask);

            // combine samples by middle-base by using degenerate code
            let mut middle_2_samples: HashMap<char, Vec<String>> = HashMap::new();
            for (i, nucl) in int_middle_base_vec.iter().enumerate() {
                // if middle-base not absent
                if *nucl != 45 {
                    for &new_nucl in &degenerate_code[nucl] {
                        middle_2_samples
                            .entry(new_nucl)
                            .or_default()
                            .push(i.to_string());
                    }
                }
            }

            for (nucl, l_indexes) in middle_2_samples.iter() {
                let str_sample_id = l_indexes.join("|");
                let value_index: u32;

                // get index to save species list (as String)
                if let Some(index) = all_indexes.get(&str_sample_id) {
                    // Key exists, use the existing index
                    value_index = *index;
                } else {
                    // Key does not exist, insert a new entry with the next index
                    value_index = all_indexes.len() as u32;
                    all_indexes.insert(str_sample_id.to_string(), value_index.clone());
                    index_map.insert(value_index.clone(), str_sample_id.to_string());
                }

                // save kmers in coloured de Bruijn graph
                let full_kmer = format!("{}{}{}", kmer_left, nucl, kmer_right);

                let encoded_kmer_1 = encode_kmer(&full_kmer[..full_kmer.len() - 1].to_string());
                let encoded_kmer_2 = encode_kmer(&full_kmer[1..].to_string());

                // uncomment to print network
                // println!("{}	{}	", &full_kmer[..full_kmer.len() - 1].to_string(), &full_kmer[1..].to_string());

                all_kmers
                    .entry(encoded_kmer_1)
                    .or_default()
                    .insert(encoded_kmer_2, value_index.clone());

                let rc_kmer = rev_compl(&full_kmer);
                let rc_encoded_kmer_1 = encode_kmer(&rc_kmer[..full_kmer.len() - 1].to_string());
                let rc_encoded_kmer_2 = encode_kmer(&rc_kmer[1..].to_string());

                // uncomment to print network
                // println!("{}	{}	", &rc_kmer[..full_kmer.len() - 1].to_string(), &rc_kmer[1..].to_string());

                all_kmers
                    .entry(rc_encoded_kmer_1)
                    .or_default()
                    .insert(rc_encoded_kmer_2, value_index);
            }
        }
    }

    // sanity check: test if the number of sample combinations is too high to be stored as u64 integers
    let max_u32_value: usize = u32::MAX as usize;
    if all_indexes.len() > max_u32_value {
        eprintln!(
            "\nError: the number of sample combinations is too high to be stored as u64 integers\n"
        );
        std::process::exit(1);
    }

    println!("     . {} nodes (includes rev-compl)", all_kmers.len());
    println!("     . {} sample combinations", all_indexes.len());
    (len_kmer, sample_names, all_kmers, index_map)
}
