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
use std::error::Error;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};

use hashbrown::HashMap;
use indicatif::ProgressIterator;
use serde::{Deserialize, Serialize};

use super::QualOpts;
use crate::io_utils::any_fastq;
use crate::ska_dict::bit_encoding::UInt;
use crate::ska_dict::SkaDict;

/// Tuple for name and fasta or paired fastq input
pub type InputFastx = (String, String, Option<String>);

/// Merged dictionary with names, and middle bases in [`Vec<u8>`] in the same order.
#[derive(Serialize, Deserialize)]
pub struct MergeSkaDict<IntT> {
    /// K-mer size
    k: usize,
    /// Whether reverse complement split k-mers were used
    rc: bool,
    /// Total number of samples supported (some may be empty)
    n_samples: usize,
    /// Sample names (some may be empty)
    names: Vec<String>,
    /// Constant bases
    ref_bases: HashMap<IntT, u8>,
    /// Dictionary of split k-mers, and middle base vectors
    variant_kmers: HashMap<IntT, Vec<(usize, u8)>>,
    added_samples: Vec<usize>,
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
            ref_bases: HashMap::default(),
            variant_kmers: HashMap::default(),
            added_samples: Vec::new(),
        }
    }

    /// Add a single [`crate::ska_dict::SkaDict`](`SkaDict`) to the merged dictionary
    ///
    /// NB: this is not a real append, and to work the input must have a unique
    /// `sample_idx < n_samples` set.
    ///
    /// # Panics
    ///
    /// If k-mer length or reverse complement do not match
    pub fn append(&mut self, other: &mut SkaDict<IntT>) {
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
        self.names[other.idx()] = other.name().to_string();

        // Vector of missing based on samples already in dict
        let mut missing_vec = Vec::new();
        for idx in self.added_samples {
            missing_vec.push((idx, b'-'));
        }

        if self.ksize() == 0 {
            self.ref_bases = other.kmers().clone();
        } else {
            // iterate over existing k-mers, and add new base or empty from
            // the new sample. Remove k-mers when processed
            for (kmer, ref_base) in self.ref_bases {
                let new_base = match other.kmers().remove(&kmer) {
                    Some(base) => {
                        if base != ref_base {
                            Some(base)
                        } else {
                            None
                        }
                    },
                    None => {
                        Some(b'-')
                    }
                };
                if let Some(base) = new_base {
                    self.variant_kmers.entry(kmer).and_modify(|vec| { vec.push((other.idx(), base))}).or_insert(vec![(other.idx(), base)]);
                }
            }
            // then iterate over remaining other and add in empty
            for (kmer, ref_base) in other.kmers() {
                self.ref_bases.insert(*kmer, *ref_base);
                self.variant_kmers.insert(*kmer, missing_vec.clone());
            }
        }
        self.added_samples.push(other.idx());
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
                swap(&mut other.added_samples, &mut self.added_samples);
                swap(&mut other.ref_bases, &mut self.ref_bases);
                swap(&mut other.variant_kmers, &mut self.variant_kmers);
            } else {
                for name_it in other.names.iter_mut().zip(self.names.iter_mut()) {
                    let (other_name, self_name) = name_it;
                    if self_name.is_empty() {
                        swap(self_name, other_name);
                    }
                }
                self.added_samples.append(&mut other.added_samples);

                let mut other_missing_vec = Vec::new();
                for idx in other.added_samples {
                    other_missing_vec.push((idx, b'-'));
                }
                let mut self_missing_vec = Vec::new();
                for idx in self.added_samples {
                    self_missing_vec.push((idx, b'-'));
                }

                let mut new_variants = HashMap::new();
                for (kmer, ref_base) in self.ref_bases {
                    // This merges the references, and deals with missing
                    match other.ref_dict().remove(&kmer) {
                        Some(other_ref) => {
                            // CASE 1 inconsistent refs must be reconciled
                            if ref_base != other_ref {
                                // For now, assuming most cases are where merging two similar sized
                                // dicts, just pick one to be ref. Better would be to pick based on MAF
                                let mut other_new = Vec::new();
                                for idx in other.added_samples {
                                    other_new.push((idx, ref_base));
                                }
                                other.variant_kmers.entry(kmer)
                                    .and_modify(|variants| {
                                        other_new.extend(variants.iter().filter(|v| v.1 != ref_base));
                                        swap(variants, &mut other_new);
                                    })
                                    .or_insert_with(|| mem::take(&mut other_new));

                            }
                            // CASE 2 else the same -- do nothing
                        },
                        // CASE 3 present in self, not in other
                        None => {
                            new_variants.insert(kmer, other_missing_vec.clone());
                        }
                    }
                }
                // CASE 4 present in other, not in self
                for (kmer, ref_base) in other.ref_dict() {
                    self.ref_bases.insert(*kmer, *ref_base);
                    new_variants.insert(*kmer, self_missing_vec.clone());
                }

                // Now deal with the variants
                // same pattern again: when in both, append
                // when only in self, do nothing
                // when only in other, move vec from other
                for (kmer, mut variants) in &other.variant_kmers {
                    self.variant_kmers.entry(*kmer)
                        .and_modify(|vec| {vec.append(&mut variants)})
                        .or_insert_with(|| mem::take(&mut variants));
                }

                // append new_variants to self.var
                self.variant_kmers.extend(new_variants);
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
        self.n_samples += other.nsamples();
        other.n_samples = self.n_samples;
        // Don't need to update names in other, they will be ignored
        self.names.extend(other.names);
        // Update indices
        for idx in &other.added_samples {
            *idx += self.n_samples;
        }
        // Also need to add to indices
        for variant in other.variant_kmers.values_mut() {
            for sample_var in variant {
                (*sample_var).0 += self.n_samples;
            }
        }

        // Then can call merge
        self.merge(other);
    }

    /// Save the split k-mer array to a `.skf` file.
    ///
    /// Serialised into [`ciborium`] format, and compressed with [`snappy`](`snap`)
    pub fn save(&self, filename: &str) -> Result<(), Box<dyn Error>> {
        let serial_file = BufWriter::new(File::create(filename)?);
        let mut compress_writer = snap::write::FrameEncoder::new(serial_file);
        ciborium::ser::into_writer(self, &mut compress_writer)?;
        Ok(())
    }

    /// Load a split k-mer array from a `.skf` file.
    pub fn load(filename: &str) -> Result<Self, Box<dyn Error>> {
        log::info!("Loading skf file");
        let ska_file = BufReader::new(File::open(filename)?);
        let decompress_reader = snap::read::FrameDecoder::new(ska_file);
        let ska_obj: Self = ciborium::de::from_reader(decompress_reader)?;
        Ok(ska_obj)
    }

    /// Write the middle bases as an alignment (FASTA).
    ///
    /// This is simply a transpose of the middle bases (which is actually
    /// recopied in memory) and written as a FASTA using [`needletail`].
    ///
    /// The writer `f` can be anything supporting the [`Write`] trait, but
    /// is typically generated by [`crate::io_utils::set_ostream`].
    ///
    /// This is used in `ska align`.
    pub fn write_fasta<W: Write>(&self, f: &mut W) -> Result<(), needletail::errors::ParseError> {
        // Do the transpose, but recopy with correct strides
        let var_t = self.variants.t();
        let mut var_t_owned = Array2::zeros(var_t.raw_dim());
        var_t_owned.assign(&var_t);

        self.names
            .iter()
            .zip(var_t_owned.outer_iter())
            .try_for_each(|it| {
                let (name, seq_u8) = it;
                write_fasta(
                    name.as_bytes(),
                    seq_u8.as_slice().expect("Array conversion error"),
                    f,
                    needletail::parser::LineEnding::Unix,
                )
            })
    }

    // TODO list
    // for first test:
    // write fasta
    // serialise IntT
    // display fn
    // map function in ska_ref
    // then:
    // debug fn
    // n_sample_kmers
    // filter
    // delete samples
    // weed
    // distance (probably no longer uses row dist)
    // (maybe) in save, recode for optimal MAF

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
    pub fn ref_dict(&self) -> &HashMap<IntT, u8> {
        &self.ref_bases
    }

    /// Total number of split k-mers
    pub fn ksize(&self) -> usize {
        self.ref_bases.len()
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
        let mut ska_dict = SkaDict::new(
            k,
            idx + offset,
            (filename, second_file.as_ref()),
            name,
            rc,
            qual,
        );
        merged_dict.append(&mut ska_dict);
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

    if any_fastq(input_files) {
        log::info!("FASTQ files filtered with: {qual}");
    } else {
        log::info!("All input files FASTA (no error filtering)");
    }

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
            let ska_dict = SkaDict::new(k, idx, (filename, second_file.as_ref()), name, rc, qual);
            merged_dict.append(&mut ska_dict);
        }
    }
    merged_dict
}
