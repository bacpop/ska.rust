use hashbrown::{HashMap, HashSet};
//use std::time::Instant;

use rayon::prelude::*;
use dashmap::DashMap;

//use crate::utils::DATA_INFO;


pub fn compact_graph(all_kmers: &mut HashMap<u128, Vec<u128>>, start_kmers: &HashSet<u128>, end_kmers: &HashSet<u128>) -> DashMap<u128, Vec<u128>> {
    
    //let data_info = DATA_INFO.get().unwrap();

    //let start = Instant::now();
        
    let compacted: DashMap<u128, Vec<u128>> = DashMap::new();
    
    // from start k-mers
    start_kmers.par_iter().for_each(|kmer| {
        if let Some(starting_kmers) = all_kmers.get(kmer) {
            for starting_kmer in starting_kmers.iter() {
                let mut current_kmer = *starting_kmer;

                let mut visited: HashSet<u128> = HashSet::new();
                let mut vec_visited: Vec<u128> = Vec::new();

                let mut walking_along_path = true;

                while walking_along_path {
                    if let Some(next_kmer_data) = all_kmers.get(&current_kmer) {
                        if next_kmer_data.len() == 1 && !visited.contains(&next_kmer_data[0]) {
                            current_kmer = next_kmer_data[0];
                            vec_visited.push(current_kmer);
                            visited.insert(current_kmer);

                            if end_kmers.contains(&current_kmer) || start_kmers.contains(&current_kmer) {
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

                let mut visited: HashSet<u128> = HashSet::new();
                let mut vec_visited: Vec<u128> = Vec::new();

                let mut walking_along_path = true;

                while walking_along_path {
                    if let Some(next_kmer_data) = all_kmers.get(&current_kmer) {
                        if next_kmer_data.len() == 1 && !visited.contains(&next_kmer_data[0]) {
                            current_kmer = next_kmer_data[0];
                            vec_visited.push(current_kmer);
                            visited.insert(current_kmer);

                            if end_kmers.contains(&current_kmer) || start_kmers.contains(&current_kmer) {
                                walking_along_path = false;
                            }
                        } else {
                            walking_along_path = false;
                        }
                    } else {
                        walking_along_path = false;
                    }
                }
                // could be "1" but for some reason I get more variant groups with k_graph
                //if vec_visited.len() > data_info.k_graph {
                if vec_visited.len() > 1 {
                    compacted.insert(*starting_kmer, vec_visited);
                }
            }
        }
    });

    //let duration = start.elapsed();
    //println!("time taken: {:?}", duration);

    //let start = Instant::now();
    
    // modify graph and compacted vector
    //let mut nb_removed = 0;
    for mut item in compacted.iter_mut() {
        let (starting_kmer, vec_visited) = item.pair_mut();
        // remove edges corresponding to compacted vector
        all_kmers.get_mut(starting_kmer).unwrap().retain(|&neighbor| neighbor != vec_visited[0]);
        for window in vec_visited[..vec_visited.len() -1].windows(2) {
            all_kmers.get_mut(&window[0]).unwrap().retain(|&neighbor| neighbor != window[1]);
            //nb_removed += 1
        }        
        
        // add new edge in place of compacted segment
        all_kmers.entry(*starting_kmer)
                 .or_default()
                 .push(vec_visited[vec_visited.len() -1]);

        // remove last element of compact vector
        vec_visited.pop();
    }
        
    //let duration = start.elapsed();
    //println!("{} edges removed", nb_removed);    
    //println!("time taken: {:?}", duration);    

    compacted
}
