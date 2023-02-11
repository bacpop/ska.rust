//! Type for combining split-kmers from multiple [`SkaDict`].
//!
//! This is an intermediate type which supports building/merging/adding, but then
//! should be converted to a [`crate::merge_ska_array::MergeSkaArray`] array for most operations.
//! Except for [`crate::ska_ref::RefSka::map()`], which implements `ska map`
//!
//! [`build_and_merge`] is the main interface.

use core::mem::swap;
use core::panic;
use std::mem;

use hashbrown::HashMap;
use indicatif::{ParallelProgressIterator, ProgressIterator};
use rayon::prelude::*;

use crate::ska_dict::SkaDict;
use crate::ska_dict::bit_encoding::RevComp;

/// Tuple for name and fasta or paired fastq input
pub type InputFastx = (String, String, Option<String>);

/// Merged dictionary with names, and middle bases in [`Vec<u8>`] in the same order.
pub struct MergeSkaDict<IntT>
where
    IntT: RevComp
{
    /// K-mer size
    k: usize,
    /// Whether reverse complement split k-mers were used
    rc: bool,
    /// Total number of samples supported (some may be empty)
    n_samples: usize,
    /// Sample names (some may be empty)
    names: Vec<String>,
    /// Dictionary of split k-mers, and middle base vectors
    split_kmers: HashMap<IntT, Vec<u8>>,
}

impl<IntT> MergeSkaDict<IntT>
where
    IntT: RevComp
{
    /// Create an empty merged dictionary, to be used with [`MergeSkaDict::merge()`]
    /// or [`MergeSkaDict::append()`].
    pub fn new(k: usize, n_samples: usize, rc: bool) -> Self {
        Self {
            k,
            rc,
            n_samples,
            names: vec!["".to_string(); n_samples],
            split_kmers: HashMap::default(),
        }
    }

    /// Directly add name and merged dictionary content
    ///
    /// Used when creating this struct from a [`crate::merge_ska_array::MergeSkaArray`]
    pub fn build_from_array<'a>(
        &'a mut self,
        names: &'a mut Vec<String>,
        split_kmers: &mut HashMap<IntT, Vec<u8>>,
    ) {
        swap(names, &mut self.names);
        swap(split_kmers, &mut self.split_kmers);
    }

    /// Add a single [`crate::ska_dict::SkaDict`](`SkaDict`) to the merged dictionary
    ///
    /// NB: this is not a real append, and to work the input must have a unique
    /// `sample_idx < n_samples` set.
    ///
    /// # Panics
    ///
    /// If k-mer length or reverse complement do not match
    pub fn append(&mut self, other: &SkaDict<IntT>) {
        if other.kmer_len() != self.k {
            panic!(
                "K-mer lengths do not match: {} {}",
                other.kmer_len(),
                self.k
            );
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

    /// Combine with another [`MergeSkaDict`] with non-overlapping samples
    ///
    /// Used when building, when individual dicts have been joined using append
    /// (CLI `ska build`)
    ///
    /// # Panics
    ///
    /// If k-mer length or reverse complement do not match
    pub fn merge<'a>(&'a mut self, other: &'a mut MergeSkaDict<IntT>) {
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
                    if self_name.is_empty() {
                        swap(self_name, other_name);
                    }
                }

                for (kmer, other_vec) in &mut other.split_kmers {
                    self.split_kmers
                        .entry(*kmer)
                        .and_modify(|self_vec| {
                            // Vectorises to VORPS (I've checked!)
                            for base_it in other_vec.iter().zip(self_vec.iter_mut()) {
                                *base_it.1 |= *base_it.0;
                            }
                        })
                        .or_insert_with(|| mem::take(other_vec));
                }
            }
        }
    }

    /// Combine with another [`MergeSkaDict`] with additional samples
    ///
    /// For joining two separate dicts (CLI 'ska merge')
    ///
    /// # Panics
    ///
    /// If k-mer length or reverse complement do not match
    pub fn extend<'a>(&'a mut self, other: &'a mut MergeSkaDict<IntT>) {
        if other.k != self.k {
            panic!("K-mer lengths do not match: {} {}", other.k, self.k);
        }
        if other.rc() != self.rc {
            panic!("Strand use inconsistent");
        }

        // Add new names in
        self.names.extend_from_slice(&other.names);
        let total_samples = self.n_samples + other.nsamples();

        // Where overlapping (and_modify) concat the base vecs
        // Where missing, add empty entries for self, and the other base vec
        for (kmer, other_vec) in &mut other.split_kmers {
            self.split_kmers
                .entry(*kmer)
                .and_modify(|self_vec| {
                    self_vec.extend_from_slice(other_vec);
                })
                .or_insert_with(|| {
                    let mut empty_samples = vec![0; self.n_samples];
                    empty_samples.extend_from_slice(other_vec);
                    empty_samples
                });
        }
        // Extend any other missing k-mers in self with empty entries for other
        for (_kmer, self_vec) in &mut self.split_kmers {
            if self_vec.len() != total_samples {
                self_vec.extend(vec![0; other.nsamples()]);
            }
        }
        self.n_samples = total_samples;
    }

    /// K-mer length of split-kmers
    pub fn kmer_len(&self) -> usize {
        self.k
    }

    /// Whether reverse complement has been used
    pub fn rc(&self) -> bool {
        self.rc
    }

    /// Sample names
    pub fn names(&self) -> &Vec<String> {
        &self.names
    }

    /// Split k-mer dictionary
    pub fn kmer_dict(&self) -> &HashMap<IntT, Vec<u8>> {
        &self.split_kmers
    }

    /// Total number of split k-mers
    pub fn ksize(&self) -> usize {
        self.split_kmers.len()
    }

    /// Number of samples
    pub fn nsamples(&self) -> usize {
        self.n_samples
    }
}

// Functions to created merged dicts from files

/// Serial `MergeSkaDict::append()` into a [`MergeSkaDict`]
fn multi_append<IntT>(
    input_dicts: &mut [SkaDict<IntT>],
    total_size: usize,
    k: usize,
    rc: bool,
) -> MergeSkaDict<IntT>
where
    IntT: RevComp
{
    let mut merged_dict = MergeSkaDict::new(k, total_size, rc);
    for ska_dict in &mut input_dicts.iter() {
        merged_dict.append(ska_dict);
    }
    merged_dict
}

/// Recursive parallel merge
///
/// Depth sets number of splits into two
/// i.e. depth 1 splits in 2, depth 2 splits in 4
fn parallel_append<IntT>(
    depth: usize,
    dict_list: &mut [SkaDict<IntT>],
    total_size: usize,
    k: usize,
    rc: bool,
) -> MergeSkaDict<IntT>
where
    IntT: RevComp
{
    let (bottom, top) = dict_list.split_at_mut(dict_list.len() / 2);
    if depth == 1 {
        let (mut bottom_merge, mut top_merge) = rayon::join(
            || multi_append(bottom, total_size, k, rc),
            || multi_append(top, total_size, k, rc),
        );
        bottom_merge.merge(&mut top_merge);
        bottom_merge
    } else {
        let (mut bottom_merge, mut top_merge) = rayon::join(
            || parallel_append(depth - 1, bottom, total_size, k, rc),
            || parallel_append(depth - 1, top, total_size, k, rc),
        );
        bottom_merge.merge(&mut top_merge);
        bottom_merge
    }
}

/// Create a [`MergeSkaDict`] from input FASTA/FASTQ files
///
/// First build [`SkaDict`] in parallel for each input.
///
/// Then merge together, in parallel if there are enough input files,
/// otherwise serially.
///
/// # Examples
/// ```
/// use ska::merge_ska_dict::{InputFastx, build_and_merge};
///
/// let input_files: [InputFastx; 2] = [("test1".to_string(),
///                                      "tests/test_files_in/test_1.fa".to_string(),
///                                      None),
///                                     ("test2".to_string(),
///                                      "tests/test_files_in/test_2.fa".to_string(),
///                                      None)];
/// let merged_dict = build_and_merge(&input_files, 17, true, 0, 0, 1);
/// ```
///
/// # Panics
///
/// If any input files are invalid
pub fn build_and_merge<IntT>(
    input_files: &[InputFastx],
    k: usize,
    rc: bool,
    min_count: u16,
    min_qual: u8,
    threads: usize,
) -> MergeSkaDict<IntT>
where
    IntT: RevComp
{
    // Build indexes
    log::info!("Building skf dicts from sequence input");
    let mut ska_dicts: Vec<SkaDict<IntT>> = Vec::new();
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
            .map(|(idx, (name, filename, second_file))| {
                SkaDict::new(
                    k,
                    idx,
                    (filename, second_file.as_ref()),
                    name,
                    rc,
                    min_count,
                    min_qual,
                )
            })
            .collect();
    } else {
        for file_it in input_files.iter().progress().enumerate() {
            let (idx, (name, filename, second_file)) = file_it;
            ska_dicts.push(SkaDict::new(
                k,
                idx,
                (filename, second_file.as_ref()),
                name,
                rc,
                min_count,
                min_qual,
            ))
        }
    }

    // Merge indexes, ensuring at least 10 samples per thread
    let mut merged_dict = MergeSkaDict::new(k, ska_dicts.len(), rc);
    let max_threads = usize::max(1, usize::min(threads, 1 + ska_dicts.len() / 10));
    let max_depth = f64::floor(f64::log2(max_threads as f64)) as usize;
    if max_depth > 0 {
        log::info!(
            "{}",
            format!(
                "Merging skf dicts in parallel using {} threads",
                1 << max_depth
            )
        );
        let total_size = ska_dicts.len();
        merged_dict = parallel_append(max_depth, &mut ska_dicts, total_size, k, rc);
    } else {
        log::info!("Merging skf dicts serially");
        for ska_dict in ska_dicts.iter_mut().progress() {
            merged_dict.append(ska_dict);
        }
    }
    merged_dict
}
