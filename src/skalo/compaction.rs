//! Functions to compact the graph in ska lo
use dashmap::DashMap;
use hashbrown::{HashMap, HashSet};
use rayon::prelude::*;

use crate::ska_dict::bit_encoding::UInt;

/// compact the assembly graph by removing all but one nodes between extremity nodes (entry or exit nodes)
/// the function outputs removed nodes in 'compacted'
pub fn compact_graph<IntT: for<'a> UInt<'a>>(
    all_kmers: &mut HashMap<IntT, Vec<IntT>>,
    start_kmers: &HashSet<IntT>,
    end_kmers: &HashSet<IntT>,
) -> DashMap<IntT, Vec<IntT>> {
    let compacted: DashMap<IntT, Vec<IntT>> = DashMap::new();

    // from start k-mers
    start_kmers.par_iter().for_each(|kmer| {
        if let Some(starting_kmers) = all_kmers.get(kmer) {
            for starting_kmer in starting_kmers.iter() {
                let mut current_kmer = *starting_kmer;

                let mut visited: HashSet<IntT> = HashSet::new();
                let mut vec_visited: Vec<IntT> = Vec::new();

                let mut walking_along_path = true;

                while walking_along_path {
                    if let Some(next_kmer_data) = all_kmers.get(&current_kmer) {
                        if next_kmer_data.len() == 1 && !visited.contains(&next_kmer_data[0]) {
                            current_kmer = next_kmer_data[0];
                            vec_visited.push(current_kmer);
                            visited.insert(current_kmer);

                            if end_kmers.contains(&current_kmer)
                                || start_kmers.contains(&current_kmer)
                            {
                                walking_along_path = false;
                            }
                        } else {
                            walking_along_path = false;
                        }
                    } else {
                        walking_along_path = false;
                    }
                }
                if vec_visited.len() > 1 {
                    compacted.insert(*starting_kmer, vec_visited);
                }
            }
        }
    });

    // from end k-mers
    end_kmers.par_iter().for_each(|kmer| {
        if let Some(starting_kmers) = all_kmers.get(kmer) {
            for starting_kmer in starting_kmers.iter() {
                let mut current_kmer = *starting_kmer;

                let mut visited: HashSet<IntT> = HashSet::new();
                let mut vec_visited: Vec<IntT> = Vec::new();

                let mut walking_along_path = true;

                while walking_along_path {
                    if let Some(next_kmer_data) = all_kmers.get(&current_kmer) {
                        if next_kmer_data.len() == 1 && !visited.contains(&next_kmer_data[0]) {
                            current_kmer = next_kmer_data[0];
                            vec_visited.push(current_kmer);
                            visited.insert(current_kmer);

                            if end_kmers.contains(&current_kmer)
                                || start_kmers.contains(&current_kmer)
                            {
                                walking_along_path = false;
                            }
                        } else {
                            walking_along_path = false;
                        }
                    } else {
                        walking_along_path = false;
                    }
                }
                if vec_visited.len() > 1 {
                    compacted.insert(*starting_kmer, vec_visited);
                }
            }
        }
    });

    // modify graph and compacted vector
    for mut item in compacted.iter_mut() {
        let (starting_kmer, vec_visited) = item.pair_mut();
        // remove edges corresponding to compacted vector
        all_kmers
            .get_mut(starting_kmer)
            .unwrap()
            .retain(|&neighbor| neighbor != vec_visited[0]);
        for window in vec_visited[..vec_visited.len() - 1].windows(2) {
            all_kmers
                .get_mut(&window[0])
                .unwrap()
                .retain(|&neighbor| neighbor != window[1]);
        }

        // add new edge in place of compacted segment
        all_kmers
            .entry(*starting_kmer)
            .or_default()
            .push(vec_visited[vec_visited.len() - 1]);

        // remove last element of compact vector
        vec_visited.pop();
    }

    compacted
}
