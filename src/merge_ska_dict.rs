// Class for split-kmers from multiple SkaDict
// In dictonary representation to support merging/adding
// Print will print out no kmers, no samples, split k-mers

use core::mem::swap;
use core::panic;
use std::fmt;
use std::mem;

use hashbrown::HashMap;
use indicatif::{ProgressIterator, ParallelProgressIterator};
use rayon::prelude::*;

use crate::ska_dict::bit_encoding::{generate_masks, decode_kmer};
use crate::ska_dict::SkaDict;
use crate::ska_ref::RefSka;

pub type InputFastx = (String, String, Option<String>);

pub struct MergeSkaDict {
    k: usize,
    rc: bool,
    n_samples: usize,
    names: Vec<String>,
    split_kmers: HashMap<u64, Vec<u8>>,
}

impl MergeSkaDict {
    pub fn kmer_len(&self) -> usize {
        self.k
    }

    pub fn rc(&self) -> bool {
        self.rc
    }

    pub fn names(&self) -> &Vec<String> {
        &self.names
    }

    pub fn kmer_dict(&self) -> &HashMap<u64, Vec<u8>> {
        &self.split_kmers
    }

    pub fn ksize(&self) -> usize {
        self.split_kmers.len()
    }

    pub fn nsamples(&self) -> usize {
        self.n_samples
    }

    pub fn new(k: usize, n_samples: usize, rc: bool) -> Self {
        let names = vec!["".to_string(); n_samples];
        let split_kmers = HashMap::default();
        return Self {
            k,
            rc,
            n_samples,
            names,
            split_kmers,
        };
    }

    pub fn build_from_array<'a>(
        &'a mut self,
        names: &'a mut Vec<String>,
        split_kmers: &mut HashMap<u64, Vec<u8>>,
    ) {
        swap(names, &mut &mut self.names);
        swap(split_kmers, &mut &mut self.split_kmers);
    }

    pub fn append(&mut self, other: &SkaDict) {
        if other.kmer_len() != self.k {
            panic!("K-mer lengths do not match: {} {}", other.kmer_len(), self.k);
        }
        if other.rc() != self.rc {
            panic!("Strand use inconsistent");
        }
        self.names[other.idx()] = other.name().clone();
        if self.ksize() == 0 {
            for (kmer, base) in other.kmers() {
                let mut base_vec: Vec<u8> = vec![0; self.n_samples];
                base_vec[other.idx()] = *base;
                self.split_kmers.insert(*kmer, base_vec);
            }
        } else {
            for (kmer, base) in other.kmers() {
                self.split_kmers
                    .entry(*kmer)
                    .and_modify(|b| {
                        b[other.idx()] = *base;
                    })
                    .or_insert_with(|| {
                        let mut new_base_vec: Vec<u8> = vec![0; self.n_samples];
                        new_base_vec[other.idx()] = *base;
                        new_base_vec
                    });
            }
        }
    }

    pub fn merge<'a>(&'a mut self, other: &'a mut MergeSkaDict) {
        if other.k != self.k {
            panic!("K-mer lengths do not match: {} {}", other.k, self.k);
        }
        if other.rc() != self.rc {
            panic!("Strand use inconsistent");
        }
        if other.ksize() > 0 {
            if self.ksize() == 0 {
                swap(&mut other.names, &mut self.names);
                swap(&mut other.split_kmers, &mut self.split_kmers);
            } else {
                for name_it in other.names.iter_mut().zip(self.names.iter_mut()) {
                    let (other_name, self_name) = name_it;
                    if self_name == "" {
                        swap(self_name, other_name);
                    }
                }

                for (kmer, other_vec) in &mut other.split_kmers {
                    self.split_kmers
                        .entry(*kmer)
                        .and_modify(|self_vec| {
                            // Vectorises to VORPS
                            for base_it in other_vec.iter().zip(self_vec.iter_mut()) {
                                *base_it.1 |= *base_it.0;
                            }
                        })
                        .or_insert_with(|| {
                            mem::take(other_vec)
                        });
                }
            }
        }
    }

    pub fn weed(&mut self, weed_ref: &RefSka) {
        let mut removed = 0;
        for kmer in weed_ref.kmer_iter() {
            if self.split_kmers.remove(&kmer).is_some() {
                removed += 1;
            }
        }
        log::info!("Removed {} of {} weed k-mers", removed, weed_ref.ksize());
    }
}

impl fmt::Display for MergeSkaDict {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "k={}\nrc={}\n{} k-mers\n{} samples\n", self.k, self.rc(), self.ksize(), self.nsamples())?;
        writeln!(f, "{:?}", self.names)
    }
}

impl fmt::Debug for MergeSkaDict {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let (lower_mask, upper_mask) = generate_masks(self.k);
        self.split_kmers.iter().try_for_each(|it| {
            let (split_kmer, vars_u8) = it;
            let mut seq_string = String::with_capacity(self.nsamples());
            for middle_base in vars_u8 {
                let base = if *middle_base == 0 {'-'} else {*middle_base as char};
                seq_string.push(base);
                seq_string.push(',');
            }
            seq_string.pop();
            let (upper, lower) = decode_kmer(self.k, *split_kmer, upper_mask, lower_mask);
            write!(f, "{}\t{}\t{}\n", upper, lower, seq_string)
        })
    }
}

// Functions to created merged dicts from files

fn multi_append(input_dicts: &mut[SkaDict], total_size: usize, k: usize, rc: bool) -> MergeSkaDict {
    let mut merged_dict = MergeSkaDict::new(k, total_size, rc);
    for ska_dict in &mut input_dicts.iter() {
        merged_dict.append(ska_dict);
    }
    return merged_dict;
}

// Recursive merge, depth sets number of splits into two i.e. depth 1 splits in 2, depth 2 splits in 4
fn parallel_append(depth: usize, dict_list: &mut [SkaDict], total_size: usize, k: usize, rc: bool) -> MergeSkaDict {
    let (bottom, top) = dict_list.split_at_mut(dict_list.len() / 2);
    if depth == 1 {
        let (mut bottom_merge, mut top_merge) =
        rayon::join(
            || multi_append(bottom, total_size, k, rc),
            || multi_append(top, total_size, k, rc)
        );
        bottom_merge.merge(&mut top_merge);
        return bottom_merge;
    } else {
        let (mut bottom_merge, mut top_merge) =
        rayon::join(
            || parallel_append(depth - 1, bottom, total_size, k, rc),
            || parallel_append(depth - 1, top, total_size, k, rc)
        );
        bottom_merge.merge(&mut top_merge);
        return bottom_merge;
    }
}

pub fn build_and_merge(
    input_files: &Vec<InputFastx>,
    k: usize,
    rc: bool,
    min_count: u16,
    min_qual: u8,
    threads: usize,
) -> MergeSkaDict {
    // Build indexes
    log::info!("Building skf dicts from sequence input");
    let mut ska_dicts: Vec<SkaDict> = Vec::new();
    ska_dicts.reserve(input_files.len());
    if threads > 1 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .unwrap();
        ska_dicts = input_files
            .par_iter()
            .progress_count(input_files.len() as u64)
            .enumerate()
            .map(|(idx, (name, filename, second_file))| SkaDict::new(k, idx, filename, second_file, name, rc, min_count, min_qual))
            .collect();
    } else {
        for file_it in input_files.iter().progress().enumerate() {
            let (idx, (name, filename, second_file)) = file_it;
            ska_dicts.push(SkaDict::new(k, idx, filename, second_file, name, rc, min_count, min_qual))
        }
    }

    // Merge indexes
    let mut merged_dict = MergeSkaDict::new(k, ska_dicts.len(), rc);
    let max_threads = usize::max(1, usize::min(threads, 1 + ska_dicts.len() / 10));
    let max_depth = f64::floor(f64::log2(max_threads as f64)) as usize;
    if max_depth > 0 {
        log::info!("{}", format!("Merging skf dicts in parallel using {} threads", 1 << max_depth));
        let total_size = ska_dicts.len();
        merged_dict = parallel_append(max_depth, &mut ska_dicts, total_size, k, rc);
    } else {
        log::info!("Merging skf dicts serially");
        for ska_dict in ska_dicts.iter_mut().progress() {
            merged_dict.append(ska_dict);
        }
    }
    return merged_dict;
}
