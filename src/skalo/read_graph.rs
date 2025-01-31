//! Path exploration in the graph
use bit_set::BitSet;
use hashbrown::{HashMap, HashSet};

use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use std::sync::{Arc, Mutex};

use crate::ska_dict::bit_encoding::UInt;
use crate::skalo::compaction::compact_graph;
use crate::skalo::process_variants::analyse_variant_groups;
use crate::skalo::utils::{Config, DataInfo, DnaSequence, VariantInfo};

/// The function traverses the graph and build variant groups (ie, set of paths starting from
/// the same entry node and ending on the same exit node) and filters them to keep only
/// paths of the most common length (unless there are only 2 paths - see indels below).
/// Indel are also identified in this function: variant groups with 2 path of diffent lengths,
/// one of which is equal or inferior to 2 * (k-1).
pub fn build_variant_groups<IntT: for<'a> UInt<'a>>(
    mut all_kmers: HashMap<IntT, Vec<IntT>>,
    start_kmers: HashSet<IntT>,
    end_kmers: HashSet<IntT>,
    kmer_2_samples: HashMap<IntT, BitSet>,
    config: &Config,
    data_info: &DataInfo,
) {
    log::info!("Compacting graph");

    let compacted = compact_graph(&mut all_kmers, &start_kmers, &end_kmers);

    log::info!("Traversing graph");

    let built_groups = Arc::new(Mutex::new(HashMap::<(IntT, IntT), Vec<VariantInfo>>::new()));

    let pool = ThreadPoolBuilder::new()
        .num_threads(config.nb_threads)
        .build()
        .unwrap();

    pool.install(|| {
        start_kmers.par_iter().for_each(|kmer| {
            let mut tmp_container: HashMap<IntT, Vec<Vec<IntT>>> = HashMap::new();

            let mut good_next: Vec<IntT> = Vec::with_capacity(2);

            for starting_kmer in all_kmers.get(kmer).unwrap().iter() {
                let mut visited = HashSet::new();
                visited.insert(*kmer);
                visited.insert(*starting_kmer);

                let mut vec_visited = vec![*kmer, *starting_kmer];

                // add compacted nodes
                if compacted.contains_key(starting_kmer) {
                    let vec_compacted = compacted.get(starting_kmer).unwrap();
                    vec_visited.extend(vec_compacted.iter());
                }

                // Initialize the stack with the starting kmer
                let mut stack = vec![PathState {
                    current_kmer: *starting_kmer,
                    visited,
                    vec_visited,
                    depth: 0,
                }];

                // Process each path in the stack
                while let Some(path_state) = stack.pop() {
                    let PathState {
                        mut current_kmer,
                        mut visited,
                        mut vec_visited,
                        depth,
                    } = path_state;

                    if depth > config.max_depth {
                        continue;
                    }

                    let mut walking_along_path = true;

                    while walking_along_path {
                        good_next.truncate(0);
                        if let Some(next_kmer_data) = all_kmers.get(&current_kmer) {
                            // add next kmers that have not yet been visited
                            for &kmer2 in next_kmer_data {
                                if !visited.contains(&kmer2) {
                                    good_next.push(kmer2);
                                }
                            }
                        }

                        match good_next.len() {
                            1 => {
                                // single path continuation
                                let next = good_next[0];
                                visited.insert(next);
                                vec_visited.push(next);
                                current_kmer = next;

                                // add compacted nodes
                                if compacted.contains_key(&next) {
                                    let vec_compacted = compacted.get(&next).unwrap();
                                    vec_visited.extend(vec_compacted.iter());
                                }

                                if end_kmers.contains(&next) {
                                    // save possible variant
                                    tmp_container
                                        .entry(next)
                                        .or_default()
                                        .push(vec_visited.clone());
                                }
                            }
                            len if len > 1 => {
                                // multiple paths -> push each onto the stack if not exit nodes
                                for next in &good_next {
                                    let mut new_visited = visited.clone();
                                    new_visited.insert(*next);

                                    let mut new_vec_visited = vec_visited.clone();
                                    new_vec_visited.push(*next);

                                    // add compacted nodes
                                    if compacted.contains_key(next) {
                                        let vec_compacted = compacted.get(next).unwrap();
                                        new_vec_visited.extend(vec_compacted.iter());
                                    }

                                    // save possible variant
                                    if end_kmers.contains(next) {
                                        tmp_container
                                            .entry(*next)
                                            .or_insert_with(Vec::new)
                                            .push(new_vec_visited.clone());
                                    }

                                    // initiate new path
                                    if walking_along_path {
                                        stack.push(PathState {
                                            current_kmer: *next,
                                            visited: new_visited,
                                            vec_visited: new_vec_visited,
                                            depth: depth + 1,
                                        });
                                    }
                                }
                                // stop current path exploration after branching
                                walking_along_path = false;
                            }
                            _ => {
                                // no further paths
                                walking_along_path = false;
                            }
                        }
                    }
                }
            }

            // save variants if at least a vector with 2+ elements for one exit k-mer
            if tmp_container.values().any(|v| v.len() > 1) {
                // prepare variant container
                let mut tmp_container_2: HashMap<(IntT, IntT), Vec<VariantInfo>> = HashMap::new();

                // check-filter-build variant groups
                for (exit_kmer, vec_variants) in tmp_container.iter() {
                    // collect second to last kmer of each variant in a hashset -> test if at least 2 (ie, the variants end on a difference)
                    let second_set: HashSet<IntT> = vec_variants.iter().map(|v| v[1]).collect();
                    let second_to_last_set: HashSet<IntT> =
                        vec_variants.iter().map(|v| v[v.len() - 2]).collect();

                    if second_set.len() > 1 && second_to_last_set.len() > 1 {
                        if let Some(most_common_length) = most_abundant_length(vec_variants) {
                            let filtered_variants: Vec<_> = if vec_variants.len() == 2 {
                                vec_variants.clone()
                            } else {
                                vec_variants
                                    .iter()
                                    .filter(|v| v.len() == most_common_length)
                                    .cloned()
                                    .collect()
                            };

                            let combined_ends: (IntT, IntT) = (*kmer, *exit_kmer);

                            // build variants 1 by 1
                            for vec_visited in filtered_variants {
                                // build sequence
                                let mut sequence = String::with_capacity(
                                    vec_visited.len() + data_info.k_graph - 1,
                                );
                                sequence
                                    .push_str(&IntT::skalo_decode_kmer(*kmer, data_info.k_graph));
                                let mut vec_snps: Vec<usize> = Vec::new();
                                for (i, next) in vec_visited.iter().enumerate() {
                                    if i != 0 {
                                        // 1st corresponds to entry k-mer
                                        sequence.push(IntT::get_last_nucl(*next));
                                    }
                                    if start_kmers.contains(next)
                                        && i <= (vec_visited.len() - data_info.k_graph)
                                    {
                                        vec_snps.push(i + data_info.k_graph);
                                    } else if end_kmers.contains(next) {
                                        vec_snps.push(i - 1);
                                    }
                                }

                                // save variant to container
                                let variant =
                                    VariantInfo::new(DnaSequence::encode(&sequence), vec_snps);

                                tmp_container_2
                                    .entry(combined_ends)
                                    .or_default()
                                    .push(variant);
                            }
                        }
                    }
                }
                // save it to main variable
                if !tmp_container_2.is_empty() {
                    let mut built_groups_locked = built_groups.lock().unwrap();
                    built_groups_locked.extend(tmp_container_2);
                }
            }
        });
    });

    let built_groups_end = built_groups.lock().unwrap();

    log::info!("{} variant groups", built_groups_end.len());

    log::info!("Identifying indels");

    // at least one of the 2 branches of an indel should have a size below or equal to this (indel and other >= (1 + 2 * data_info.k_graph))
    let min_indel = 2 * data_info.k_graph;

    // separate indels from the other variants
    let mut final_groups: HashMap<(IntT, IntT), Vec<VariantInfo>> = HashMap::new();
    let mut final_indels: HashMap<(IntT, IntT), Vec<VariantInfo>> = HashMap::new();

    for (extremities_combined, vec_variant) in built_groups_end.iter() {
        // test if variant is an indel
        if vec_variant.len() < 2 {
            continue;
        } else if vec_variant.len() == 2
            && vec_variant[0].sequence.len() != vec_variant[1].sequence.len()
        {
            let mut is_indel = false;
            for variant in vec_variant {
                if variant.sequence.len() <= min_indel {
                    is_indel = true;
                }
            }
            if is_indel {
                final_indels.insert(*extremities_combined, vec_variant.clone());
            }
        } else {
            final_groups.insert(*extremities_combined, vec_variant.clone());
        }
    }

    // infer variants
    analyse_variant_groups(
        final_groups,
        final_indels,
        kmer_2_samples,
        config,
        data_info,
    );
}

// find the most abundant length in a vector of variants
fn most_abundant_length<IntT: for<'a> UInt<'a>>(vec_variants: &[Vec<IntT>]) -> Option<usize> {
    let mut length_counts = std::collections::HashMap::new();

    // count the frequency of each length
    for variant in vec_variants {
        *length_counts.entry(variant.len()).or_insert(0) += 1;
    }

    // find the length with the maximum count
    length_counts
        .into_iter()
        .max_by_key(|&(_, count)| count)
        .map(|(length, _)| length)
}

/// structure to hold state for each path in the stack
pub struct PathState<T> {
    current_kmer: T,
    visited: HashSet<T>,
    vec_visited: Vec<T>,
    depth: usize,
}
