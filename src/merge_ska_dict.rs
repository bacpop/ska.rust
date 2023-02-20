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
use indicatif::ProgressIterator;

use super::QualOpts;
use crate::ska_dict::bit_encoding::UInt;
use crate::ska_dict::SkaDict;
use crate::ska_dict::count_min_filter::{CountMin, KmerFilter, BlockedBloom};

/// Tuple for name and fasta or paired fastq input
pub type InputFastx = (String, String, Option<String>);

/// Merged dictionary with names, and middle bases in [`Vec<u8>`] in the same order.
pub struct MergeSkaDict<IntT> {
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
    IntT: for<'a> UInt<'a>,
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
    pub fn append<FiltT>(&mut self, other: &SkaDict<IntT, FiltT>)
    where
        FiltT: KmerFilter
    {
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
    input_files: &[InputFastx],
    offset: usize,
    total_size: usize,
    k: usize,
    rc: bool,
    qual: &QualOpts,
) -> MergeSkaDict<IntT>
where
    IntT: for<'a> UInt<'a>,
{
    let mut merged_dict = MergeSkaDict::new(k, total_size, rc);
    for (idx, (name, filename, second_file)) in input_files.iter().enumerate() {
        if qual.bloom() {
            let filter = BlockedBloom::empty(qual.min_count);
            let ska_dict = SkaDict::new(
                k,
                idx + offset,
                (filename, second_file.as_ref()),
                name,
                rc,
                qual,
                &filter,
            );
            merged_dict.append(&ska_dict);
        } else {
            let filter = CountMin::empty(qual.min_count);
            let ska_dict = SkaDict::new(
                k,
                idx + offset,
                (filename, second_file.as_ref()),
                name,
                rc,
                qual,
                &filter,
            );
            merged_dict.append(&ska_dict);
        }

    }
    merged_dict
}

/// Recursive parallel merge
///
/// Depth sets number of splits into two
/// i.e. depth 1 splits in 2, depth 2 splits in 4
fn parallel_append<IntT>(
    depth: usize,
    offset: usize,
    file_list: &[InputFastx],
    total_size: usize,
    k: usize,
    rc: bool,
    qual: &QualOpts,
) -> MergeSkaDict<IntT>
where
    IntT: for<'a> UInt<'a>,
{
    let split_point = file_list.len() / 2;
    let (bottom, top) = file_list.split_at(split_point);
    if depth == 1 {
        let (mut bottom_merge, mut top_merge) = rayon::join(
            || multi_append(bottom, offset, total_size, k, rc, qual),
            || multi_append(top, offset + split_point, total_size, k, rc, qual),
        );
        bottom_merge.merge(&mut top_merge);
        bottom_merge
    } else {
        let (mut bottom_merge, mut top_merge) = rayon::join(
            || parallel_append(depth - 1, offset, bottom, total_size, k, rc, qual),
            || {
                parallel_append(
                    depth - 1,
                    offset + split_point,
                    top,
                    total_size,
                    k,
                    rc,
                    qual,
                )
            },
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
/// use ska::{QualOpts, QualFilter};
///
/// let quality = QualOpts {min_count: 1, min_qual: 0, qual_filter: QualFilter::NoFilter};
/// let input_files: [InputFastx; 2] = [("test1".to_string(),
///                                      "tests/test_files_in/test_1.fa".to_string(),
///                                      None),
///                                     ("test2".to_string(),
///                                      "tests/test_files_in/test_2.fa".to_string(),
///                                      None)];
/// let merged_dict = build_and_merge::<u64>(&input_files, 17, true, &quality, 1);
/// ```
///
/// # Panics
///
/// If any input files are invalid
pub fn build_and_merge<IntT>(
    input_files: &[InputFastx],
    k: usize,
    rc: bool,
    qual: &QualOpts,
    threads: usize,
) -> MergeSkaDict<IntT>
where
    IntT: for<'a> UInt<'a>,
{
    // Build indexes
    log::info!("Building skf dicts from sequence input");
    log::info!("If FASTQ input: {qual}");
    if threads > 1 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .unwrap();
    }

    // Merge indexes, ensuring at least 10 samples per thread
    let total_size = input_files.len();
    let mut merged_dict = MergeSkaDict::new(k, total_size, rc);
    let max_threads = usize::max(1, usize::min(threads, 1 + total_size / 10));
    let max_depth = f64::floor(f64::log2(max_threads as f64)) as usize;
    if max_depth > 0 {
        log::info!(
            "{}",
            format!(
                "Build and merge skf dicts in parallel using {} threads",
                1 << max_depth
            )
        );
        merged_dict = parallel_append(max_depth, 0, input_files, total_size, k, rc, qual);
    } else {
        log::info!("Build and merge serially");
        for (idx, (name, filename, second_file)) in input_files.iter().progress().enumerate() {
            if qual.bloom() {
                let filter = BlockedBloom::empty(qual.min_count);
                let ska_dict = SkaDict::new(k, idx, (filename, second_file.as_ref()), name, rc, qual, &filter);
                merged_dict.append(&ska_dict);
            } else {
                let filter = CountMin::empty(qual.min_count);
                let ska_dict = SkaDict::new(k, idx, (filename, second_file.as_ref()), name, rc, qual, &filter);
                merged_dict.append(&ska_dict);
            }
        }
    }
    merged_dict
}
