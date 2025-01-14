use hashbrown::{HashMap, HashSet};
use std::str::FromStr;

use crate::skalo_utils::{decode_kmer, get_last_nucleotide, rev_compl_u128, VariantInfo};

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
            let mut sequence = decode_kmer(kmer.clone(), len_kmer_graph);
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

fn jaccard_similarity(set1: &HashSet<&str>, set2: &HashSet<&str>) -> f32 {
    // returns the Jaccard similarity value between 2 sets
    let intersection_size = set1.intersection(set2).count() as f32;
    let union_size = (set1.len() as f32 + set2.len() as f32 - intersection_size) as f32;
    intersection_size / union_size
}

fn count_occurrences(vec: &Vec<u32>) -> HashMap<u32, i32> {
    let mut count_map = HashMap::new();
    for num in vec {
        *count_map.entry(*num).or_insert(0) += 1;
    }
    count_map
}

fn convert_combined(input: &str, k: usize) -> Result<String, String> {
    let parts: Vec<&str> = input.split('@').collect();

    // parse the u128 values from the string parts
    let u128_1 = u128::from_str(parts[0]).map_err(|e| e.to_string())?;
    let u128_2 = u128::from_str(parts[1]).map_err(|e| e.to_string())?;

    // compute the reverse complements for both u128 values using the same k
    let rev_compl_1 = rev_compl_u128(u128_1, k);
    let rev_compl_2 = rev_compl_u128(u128_2, k);

    Ok(format!("{}@{}", rev_compl_2, rev_compl_1))
}
