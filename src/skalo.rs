use crate::MergeSkaArray;
use hashbrown::{HashMap, HashSet};
use std::fs::File;
use std::io::Write;

use crate::ska_dict::bit_encoding::decode_kmer;
use crate::skalo_utils::{
    calculate_ratio_missing, check_homopolymer, collect_middle_bases, compare_samples,
    convert_combined, count_occurrences, decode_kmer_skalo, encode_kmer, get_last_nucleotide,
    jaccard_similarity, rev_compl, test_multiple_positions, VariantInfo,
};

/// Reads .skf files for skalo
pub fn ska_array_to_dbgraph(
    ska_array: MergeSkaArray<u128>,
) -> (
    usize,
    Vec<String>,
    HashMap<u128, HashMap<u128, u32>>,
    HashMap<u32, String>,
) {
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

pub fn identify_good_kmers(
    len_kmer: usize,
    all_kmers: &HashMap<u128, HashMap<u128, u32>>,
    index_map: &HashMap<u32, String>,
) -> (HashSet<u128>, HashSet<u128>) {
    println!(" # identify bubble extremities");

    let len_kmer_graph = len_kmer - 1;

    let mut start_kmers: HashSet<u128> = HashSet::new();
    let mut end_kmers: HashSet<u128> = HashSet::new();

    for (kmer, next_kmers_map) in all_kmers.iter() {
        if next_kmers_map.len() > 1 {
            let all_next_kmer: Vec<&u128> = next_kmers_map.keys().collect();

            'i_loop: for (i, &next_kmer) in all_next_kmer.iter().enumerate() {
                for &next_kmer2 in all_next_kmer.iter().skip(i + 1) {
                    let samples_1 = &index_map[&next_kmers_map[next_kmer]];
                    let samples_2 = &index_map[&next_kmers_map[next_kmer2]];

                    if compare_samples(samples_1, samples_2) {
                        start_kmers.insert(kmer.clone());

                        let dna = decode_kmer_skalo(kmer.clone(), len_kmer_graph);
                        let rc = rev_compl(&dna);

                        //uncomment to print network
                        //println!("{}	{}	red", &dna, &dna);

                        end_kmers.insert(encode_kmer(&rc));

                        //uncomment to print network
                        //println!("{}	{}	red", &rc, &rc);
                        break 'i_loop;
                    }
                }
            }
        }
    }

    // exit program if no extremity found (e.g. cases of weeded skf files)
    if start_kmers.is_empty() {
        eprintln!("\n      Error: there is no entry node in this graph, hence no variant.\n");
        std::process::exit(1);
    }

    println!("     . {} entry nodes", start_kmers.len());
    (start_kmers, end_kmers)
}

pub fn build_sequences<'a>(
    len_kmer: usize,
    all_kmers: &HashMap<u128, HashMap<u128, u32>>,
    start_kmers: &HashSet<u128>,
    end_kmers: &HashSet<u128>,
    index_map: &'a HashMap<u32, String>,
) -> HashMap<String, Vec<VariantInfo<'a>>> {
    println!(" # explore graph");

    let len_kmer_graph = len_kmer - 1;
    let len_snp = (2 * len_kmer_graph) + 1;

    let mut built_groups: HashMap<String, Vec<VariantInfo>> = HashMap::new();

    // scan for variants from entry nodes
    for kmer in start_kmers {
        // get k-mers to test
        let mut kmers_to_test = Vec::<&u128>::new();
        for (next_kmer, _) in all_kmers.get(kmer).unwrap().iter() {
            kmers_to_test.push(next_kmer);
        }
        kmers_to_test.sort();

        // consider all nucleotides and build sequence for each of them
        for next_kmer in kmers_to_test {
            // create new sequence by adding last nucleotide of the next kmer to the main kmer
            let mut sequence = decode_kmer_skalo(kmer.clone(), len_kmer_graph);
            let next_nucl = get_last_nucleotide(next_kmer.clone());
            sequence += &next_nucl.to_string();

            // create sample variables to be used during the walk
            let tmp_index_samples = all_kmers[kmer][next_kmer];
            let sample_set_ref: HashSet<&str> = index_map
                .get(&tmp_index_samples)
                .unwrap()
                .split('|')
                .collect();

            let mut d_nb_samples: HashMap<&str, i32> = HashMap::new();
            for sample in sample_set_ref.iter() {
                *d_nb_samples.entry(sample).or_insert(0) += 1;
            }

            // create set of visited nodes
            let mut visited: HashSet<u128> = HashSet::new();
            visited.insert(kmer.clone());
            visited.insert(next_kmer.clone());

            let mut previous_kmer: u128 = next_kmer.clone();
            let mut walking_along_path = true;

            let mut tmp_l_index = Vec::new();
            tmp_l_index.push(tmp_index_samples);

            while walking_along_path {
                // compare samples of all next kmers with s_ref_samples
                let mut good_next = Vec::new();

                if let Some(next_kmer_data) = all_kmers.get(&previous_kmer) {
                    // get among next kmer(s) those that were not already visited
                    let mut kmer2_to_test = Vec::new();
                    for kmer2 in next_kmer_data.keys() {
                        // we ignore already visited kmers to avoid endless loops
                        if !visited.contains(kmer2) {
                            kmer2_to_test.push(kmer2.clone());
                        }
                    }

                    // if only 1 kmer -> save it; if not -> compute Jaccard similarity
                    if kmer2_to_test.len() == 1 {
                        good_next.extend_from_slice(&kmer2_to_test);
                    } else {
                        // /*  uncomment here and below to remove the path selection based on Jaccard similarity
                        let mut max_sim = 0.0;
                        for kmer2 in kmer2_to_test {
                            let index_samples = next_kmer_data.get(&kmer2).unwrap();
                            let s_samples: HashSet<&str> =
                                index_map.get(index_samples).unwrap().split('|').collect();
                            let sim = jaccard_similarity(&sample_set_ref, &s_samples);
                            if sim > max_sim {
                                good_next.clear();
                                good_next.push(kmer2.clone());
                                max_sim = sim;
                            } else if sim == max_sim {
                                good_next.push(kmer2.clone());
                            }
                        }
                        // */
                    }
                }

                // case only 1 next kmer
                if good_next.len() == 1 {
                    // update sequence
                    sequence.push(get_last_nucleotide(good_next[0]));

                    // update visited nodes
                    visited.insert(good_next[0]);

                    // update visited sample combination indexes
                    tmp_l_index.push(all_kmers[&previous_kmer][&good_next[0]]);

                    // update previous_kmer
                    previous_kmer = good_next[0].clone();

                    // save sequence if the kmer was an end kmer
                    if end_kmers.contains(&good_next[0]) {
                        // get consensus limit to rebuild s_ref_samples
                        let limit_consensus = (visited.len() as f32 * 0.5) as i32;

                        // build majrule_samples with majority rule consensus
                        for (tmp_index, tmp_count) in count_occurrences(&tmp_l_index) {
                            let tmp_samples: HashSet<&str> =
                                index_map.get(&tmp_index).unwrap().split('|').collect();
                            for sample in tmp_samples.iter() {
                                d_nb_samples
                                    .entry(sample)
                                    .and_modify(|count| *count += tmp_count)
                                    .or_insert(tmp_count);
                            }
                        }
                        tmp_l_index.clear();

                        let mut majrule_samples: HashSet<&str> = HashSet::new();
                        for (sample, nb) in &d_nb_samples {
                            if nb >= &limit_consensus {
                                majrule_samples.insert(&sample);
                            }
                        }

                        // save variant to variant group
                        let combined_ends = format!("{}@{}", kmer, good_next[0]);

                        let variant = VariantInfo::new(
                            //kmer.clone(),
                            //good_next[0].clone(),
                            sequence.clone(),
                            false,
                            visited.clone(),
                            d_nb_samples.clone(),
                            majrule_samples,
                        );

                        built_groups
                            .entry(combined_ends.clone())
                            .or_insert_with(Vec::new)
                            .push(variant);

                        // stop looking if another branch already ended on this exit kmer
                        if built_groups[&combined_ends].len() > 1 {
                            walking_along_path = false;
                        }
                    }
                // case where no next kmer or several possible next kmers
                } else {
                    walking_along_path = false;
                }
            }
        }
    }

    // process variant groups
    let mut final_groups: HashMap<String, Vec<VariantInfo>> = HashMap::new();
    let mut nb_snps = 0;

    for (extremities_combined, vec_variant) in built_groups.iter_mut() {
        //only consider groups with 2+ variants
        if vec_variant.len() > 1 {
            // remove duplicated groups (reverse-complement) by selecting the group with highest number of maj-rule samples (-> stability of output)
            let mut to_save = false;
            let mut to_replace = false;

            let rc_extremities = convert_combined(extremities_combined, len_kmer_graph).unwrap();

            if final_groups.contains_key(&rc_extremities) {
                let nb1_majrule_samples: usize = vec_variant
                    .iter()
                    .map(|variant| variant.maj_samples.len())
                    .sum();
                let vec_variant2 = final_groups.get(&rc_extremities).unwrap();
                let nb2_majrule_samples: usize = vec_variant2
                    .iter()
                    .map(|variant| variant.maj_samples.len())
                    .sum();

                if nb1_majrule_samples > nb2_majrule_samples {
                    to_replace = true;
                    // delete rev-compl group
                    final_groups.remove(&rc_extremities);
                }
            } else {
                to_save = true;
            }

            if to_save || to_replace {
                // check if variant group is a SNP
                let seq_lengths: Vec<usize> = vec_variant
                    .iter()
                    .map(|variant| variant.sequence.len())
                    .collect();
                if seq_lengths.iter().all(|&len| len == len_snp) {
                    // update all variants
                    for variant in vec_variant.iter_mut() {
                        variant.is_snp = true;
                    }
                    // update number of SNPS if not replacement
                    if to_save {
                        nb_snps += 1;
                    }
                }

                // sort vector of variants by their number of samples (decreasing order)
                vec_variant.sort_by(|a, b| b.maj_samples.len().cmp(&a.maj_samples.len()));

                // finally save the variant group
                final_groups.insert(extremities_combined.clone(), vec_variant.clone());
            }
        }
    }

    println!("     . {} snps", nb_snps);
    println!("     . {} indels/complex", final_groups.len() - nb_snps);
    final_groups
}

pub fn filter_output_sequences(
    variant_groups: HashMap<String, Vec<VariantInfo>>,
    len_kmer: usize,
    all_samples: Vec<String>,
    n_homopolymer: Option<u32>,
    max_missing: f32,
    output_name: &str,
    input_name: &str,
) {
    println!(" # filter sequences");

    let len_kmer_graph = len_kmer - 1;

    // prepare output files
    let mut out_1 =
        File::create(format!("{}_seq_groups.fas", output_name)).expect("Failed to create file");
    let mut out_2 =
        File::create(format!("{}_variants.tsv", output_name)).expect("Failed to create file");

    // write header of TSV variant file
    out_2
        .write_all(format!("#skalo version {}\n", env!("CARGO_PKG_VERSION")).as_bytes())
        .expect("Failed to write to file");
    out_2
        .write_all(
            format!(
                "#parameters: m={}{}\n",
                max_missing,
                match n_homopolymer {
                    Some(n) => format!(" n={}", n),
                    None => String::new(),
                }
            )
            .as_bytes(),
        )
        .expect("Failed to write to file");
    out_2
        .write_all(format!("#input file: {} (k={})\n", input_name, len_kmer).as_bytes())
        .expect("Failed to write to file");

    let index_name: Vec<String> = all_samples
        .iter()
        .enumerate()
        .map(|(index, value)| format!("{}:{}", index, value))
        .collect();
    let formatted_string = index_name.join(", ");
    out_2
        .write_all(format!("#samples: {}\n", formatted_string).as_bytes())
        .expect("Failed to write to file");
    out_2.write_all(b"#pos_ali\tnb_states\ttype\tunclear_insert\tfirst_kmer\tvariants\tlast_kmer\tratio_missing\tsamples\n")
        .expect("Failed to write to file");

    // Prepare variables for binary alignment
    let mut binary_seq: HashMap<String, Vec<String>> = HashMap::new();
    let mut d_samples: HashMap<String, String> = HashMap::new();
    for (i, sample) in all_samples.iter().enumerate() {
        binary_seq.insert(sample.clone(), Vec::new());
        d_samples.insert(i.to_string(), sample.clone()); // map sample ID to sample full name
    }

    let mut nb_missing = 0;
    let mut nb_homopolymer = 0;

    // variant groups 1 by 1
    let mut position = 0;
    for vec_variants in variant_groups.values() {
        if vec_variants[0].is_snp {

            // do something with SNPs
        } else {
            // check the ratio of missing samples
            let ratio_missing = calculate_ratio_missing(all_samples.len(), vec_variants.clone());

            if ratio_missing > max_missing {
                nb_missing += 1;
            } else {
                // CHECK IF UNCLEAR INSERT
                // extract 1st kmer using 1st sequence
                let first_kmer = vec_variants[0].sequence[..len_kmer_graph].to_string();

                // collect 'middle-bases' and last k-mer (forward and reverse-complement)
                let (l_middle_bases, last_kmer) =
                    collect_middle_bases(vec_variants.clone(), len_kmer_graph, false);
                let (rc_l_middle_bases, rc_last_kmer) =
                    collect_middle_bases(vec_variants.clone(), len_kmer_graph, true);

                // test if multiple positions exist for the last k-mer
                let multiple_positions =
                    test_multiple_positions(l_middle_bases.clone(), last_kmer.clone());
                let multiple_positions_rc =
                    test_multiple_positions(rc_l_middle_bases, rc_last_kmer);

                let mut unclear_insert = "no";
                if multiple_positions || multiple_positions_rc {
                    unclear_insert = "yes";
                }

                // CHECK IK HOMOPOLYMER
                // check if indel in homopolymer
                let mut is_homopolymer = false;
                if let Some(mut n_max) = n_homopolymer {
                    // adjust n_homopolymer if higher than k-mer length to avoid the program crashing
                    if n_max > len_kmer_graph.try_into().unwrap() {
                        n_max = len_kmer_graph as u32;
                        println!("     . n_homopolymer was reduced to fit the k-mer length");
                    }
                    // test presence homopolymer
                    is_homopolymer =
                        check_homopolymer(n_max, first_kmer.clone(), l_middle_bases.clone());
                }

                if is_homopolymer {
                    nb_homopolymer += 1;
                } else {
                    // Output seq_groups
                    for variant in vec_variants.iter() {
                        let str_samples = variant
                            .maj_samples
                            .iter()
                            .copied()
                            .collect::<Vec<_>>()
                            .join("|");
                        out_1
                            .write_all(
                                format!(">{}_{}\n{}\n", position, str_samples, variant.sequence)
                                    .as_bytes(),
                            )
                            .expect("Failed to write to file");
                    }

                    // update binary alignment ('-' = missing data; '?' = both states (= 'unknown'))
                    // only consider the 2 most frequent variants in cases of 3+
                    let mut sample_done: HashSet<String> = HashSet::new();
                    // let mut state = 0;
                    for (state, variant) in vec_variants.iter().take(2).enumerate() {
                        // collect sample names for this sequence and update its vector in binary_seq
                        for sample_id in &variant.maj_samples {
                            let full_name = d_samples.get(sample_id.to_owned()).unwrap();
                            sample_done.insert(full_name.clone());
                            let sample_vec = binary_seq.get_mut(full_name).unwrap();

                            if sample_vec.len() <= position {
                                sample_vec.resize(position + 1, state.to_string());
                            } else {
                                sample_vec[position] = "?".to_string();
                            }
                        }
                        // update number of character states
                        // state += 1;
                    }

                    // update binary alignment with missing samples
                    for sample in &all_samples {
                        if !sample_done.contains(sample) {
                            binary_seq
                                .entry(sample.to_string())
                                .and_modify(|vec| vec.push("-".to_string()));
                        }
                    }

                    // get type of variant
                    let mut type_variant = "complex";
                    for xx in &l_middle_bases {
                        if xx == "." {
                            type_variant = "_indel_";
                        }
                    }

                    // get lists of samples
                    let mut l_samples = Vec::new();
                    for variant in vec_variants.iter() {
                        // sort list of samples by equivalent integers
                        let mut vect_idx = variant.maj_samples.iter().copied().collect::<Vec<_>>();
                        vect_idx.sort_by(|a, b| {
                            let a_int = a.parse::<i32>().unwrap_or(i32::MAX);
                            let b_int = b.parse::<i32>().unwrap_or(i32::MAX);
                            a_int.cmp(&b_int)
                        });
                        // save it
                        l_samples.push(vect_idx.join(","));
                    }

                    // round ratio_missing to the 2nd decimal
                    let rounded_missing = (ratio_missing * 100.0).round() / 100.0;

                    // save variants information
                    out_2
                        .write_all(
                            format!(
                                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                                position,
                                l_middle_bases.len(),
                                type_variant,
                                unclear_insert,
                                first_kmer,
                                l_middle_bases.join(" / "),
                                last_kmer,
                                rounded_missing,
                                l_samples.join(" / ")
                            )
                            .as_bytes(),
                        )
                        .expect("Failed to write to file");

                    // Update position
                    position += 1;
                }
            }
        }
    }

    // Close output files
    out_1.flush().expect("Failed to flush file");
    out_2.flush().expect("Failed to flush file");

    // Output alignment
    let mut out_3 =
        File::create(format!("{}_binary_ali.fas", output_name)).expect("Failed to create file");
    for (sample, l_seq) in &binary_seq {
        out_3
            .write_all(format!(">{}\n{}\n", sample, l_seq.join("")).as_bytes())
            .expect("Failed to write to file");
    }
    out_3.flush().expect("Failed to flush file");

    // final message
    println!("     . {} removed because of missing data", nb_missing);
    if n_homopolymer.is_some() {
        println!("     . {} removed because in homopolymers", nb_homopolymer);
    }
    println!("     . {} FINAL variant groups", position);
    println!("done.");
}
