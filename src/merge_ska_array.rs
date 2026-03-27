//! Main type for working with multiple samples.
//!
//! In fixed size array representation to support:
//!  - filter
//!  - save/load
//!  - distances
//!  - mapping
//!  - sample deletion
//!  - printing out alignment
//!
//! Can be converted to/from [`MergeSkaDict`], which is also how to create
//! from FASTA/FASTQ input.
//!
//! [`fmt::Display`] and [`fmt::Debug`] implemented to support `ska nk`.

use core::panic;
use std::error::Error;
use std::fmt;
use std::fs::File;
#[cfg(target_arch = "wasm32")]
use std::io::{BufReader, BufWriter};
#[cfg(not(target_arch = "wasm32"))]
use std::io::{BufReader, BufWriter, Write};
use std::mem;

use hashbrown::{HashMap, HashSet};
use indicatif::ParallelProgressIterator;
use ndarray::parallel::prelude::*;
use ndarray::{Array2, ArrayView, Axis, Dim};
#[cfg(not(target_arch = "wasm32"))]
use needletail::parser::write_fasta;
use rayon::iter::ParallelIterator;
use serde::{Deserialize, Serialize};

#[cfg(target_arch = "wasm32")]
use crate::logw;
use crate::merge_ska_dict::MergeSkaDict;
use crate::ska_dict::bit_encoding::{base_to_prob, decode_kmer, is_ambiguous, UInt};
use crate::ska_ref::RefSka;

#[cfg(not(target_arch = "wasm32"))]
use crate::cli::FilterType;

/// Contains information on ska distance results
pub struct VariantDist {
    /// Number of SNPs different (may be non-integer due to ambig bases)
    distance: f64,
    /// Proportion of k-mers mismatching
    mismatch_prop: f64,
    /// Absolute count of matching sites
    match_count: usize,
    /// Absolute count of mismatching sites
    mismatch_count: usize,
}

/// Prints distance, mismatch proportion, count matching, count mismatching (tab separated)
impl fmt::Display for VariantDist {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{:.2}\t{:.5}\t{}\t{}",
            self.distance, self.mismatch_prop, self.match_count, self.mismatch_count,
        )
    }
}

/// Array representation of split k-mers from multiple samples.
///
/// Supports most modification and input/output.
///
/// # Examples
/// ```
/// use ska::merge_ska_array::MergeSkaArray;
/// use ska::io_utils::set_ostream;
/// use ska::ska_ref::RefSka;
/// use ska::cli::FilterType;
///
/// // Load an array from file
/// let mut ska_array = MergeSkaArray::<u64>::load(&"tests/test_files_in/merge.skf").expect("Could not open array");
///
/// // Write alignment as FASTA on stdout
/// let mut alignment_file = set_ostream(&None);
/// ska_array.write_fasta(&mut alignment_file);
///
/// // Remove constant sites and save
/// let min_count = 1;                          // no filtering by minor allele frequency
/// let filter_ambig_as_missing = false;        // allow ambiguous bases when counting allele frequency
/// let filter = FilterType::NoAmbigOrConst;    // remove sites with no minor allele
/// let mask_ambiguous = false;                 // leave ambiguous sites as IUPAC codes
/// let ignore_const_gaps = false;              // keep sites with only '-' as variants
/// let update_counts = true;                   // keep counts updated, as saving
/// ska_array.filter(min_count, filter_ambig_as_missing, &filter, mask_ambiguous, ignore_const_gaps, update_counts);
/// ska_array.save(&"no_const_sites.skf");
///
/// // Create an iterators
/// let mut kmer_iter = ska_array.iter();
/// let (kmer, middle_base_vec) = kmer_iter.next().unwrap();
///
/// // Delete a sample
/// ska_array.delete_samples(&[&"test_1"]);
///
/// // Delete k-mers
/// let mask_repeats = false;
/// let ska_weed = RefSka::new(ska_array.kmer_len(), &"tests/test_files_in/weed.fa", ska_array.rc(), mask_repeats, mask_ambiguous);
/// let reverse = false;
/// ska_array.weed(&ska_weed, reverse);
/// ```
#[derive(Serialize, Deserialize)]
pub struct MergeSkaArray<IntT> {
    /// K-mer size
    k: usize,
    /// Whether reverse complement split k-mers were used
    rc: bool,
    /// Sample names
    names: Vec<String>,
    /// List of split k-mers
    split_kmers: Vec<IntT>,
    /// Array of middle bases, rows same order as split k-mers, columns same order as names
    variants: Array2<u8>,
    /// Count of non-missing bases for each split k-mer
    variant_count: Vec<usize>,
    /// SKA version
    ska_version: String,
    /// k-mer bits
    k_bits: u32,
}

impl<IntT> MergeSkaArray<IntT>
where
    IntT: for<'a> UInt<'a>,
{
    /// Update `variant_count` after changing `variants`.
    ///
    /// Recalculates counts, and removes any totally empty rows.
    ///
    /// # Arguments
    ///
    /// - `filter_ambig_as_missing` -- any non-ACGTU base counts as missing.
    fn update_counts(&mut self, filter_ambig_as_missing: bool) {
        log::info!("Updating variant counts");
        let mut new_counts = Vec::with_capacity(self.variant_count.len());
        let mut new_sk = Vec::with_capacity(self.split_kmers.len());

        let mut empty: usize = 0;
        let mut new_variants = Array2::zeros((0, self.names.len()));
        for (var_row, sk) in self.variants.outer_iter().zip(self.split_kmers.iter()) {
            let count = var_row
                .iter()
                .filter(|b| **b != b'-' && (!filter_ambig_as_missing || !is_ambiguous(**b)))
                .count();
            if count > 0 {
                new_counts.push(count);
                new_sk.push(*sk);
                new_variants.push_row(var_row).unwrap();
            } else {
                empty += 1;
            }
        }
        log::info!("Removed {empty} empty variants");
        self.split_kmers = new_sk;
        self.variants = new_variants;
        self.variant_count = new_counts;
    }

    /// Convert a dynamic [`MergeSkaDict`] to static array representation.
    pub fn new(dynamic: &MergeSkaDict<IntT>) -> Self {
        let mut variants = Array2::zeros((0, dynamic.nsamples()));
        let mut split_kmers: Vec<IntT> = Vec::with_capacity(dynamic.ksize());
        let mut variant_count: Vec<usize> = Vec::new();
        for (kmer, bases) in dynamic.kmer_dict() {
            split_kmers.push(*kmer);
            variant_count.push(bases.iter().filter(|b| **b != 0 && **b != b'-').count());
            variants.push_row(ArrayView::from(bases)).unwrap();
        }
        variants.mapv_inplace(|b| u8::max(b, b'-')); // turns zeros to missing
        Self {
            k: dynamic.kmer_len(),
            rc: dynamic.rc(),
            names: dynamic.names().clone(),
            split_kmers,
            variants,
            variant_count,
            ska_version: env!("CARGO_PKG_VERSION").to_string(),
            k_bits: IntT::n_bits(),
        }
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
        let ska_file = BufReader::new(File::open(filename)?);
        let decompress_reader = snap::read::FrameDecoder::new(ska_file);
        let ska_obj: Self = ciborium::de::from_reader(decompress_reader)?;
        Ok(ska_obj)
    }

    /// Convert to a dictionary representation [`MergeSkaDict`].
    ///
    /// Necessary if adding more samples.
    pub fn to_dict(&self) -> MergeSkaDict<IntT> {
        let n_samples = self.names.len();
        let mut names = self.names.clone();
        let mut split_kmers: HashMap<IntT, Vec<u8>> = HashMap::new();
        split_kmers.reserve(self.variants.nrows());
        for row_it in self.variants.outer_iter().zip(self.split_kmers.iter()) {
            let (row_vec, kmer) = row_it;
            split_kmers.insert(*kmer, row_vec.to_vec());
        }
        let mut dict = MergeSkaDict::new(self.k, n_samples, self.rc);
        dict.build_from_array(&mut names, &mut split_kmers);
        dict
    }

    /// Delete a list of named samples.
    ///
    /// Also updates counts and removes any completely missing k-mers.
    ///
    /// # Panics
    ///
    /// - If any input sample name is not in the array.
    /// - If no samples or all samples are being removed.
    pub fn delete_samples(&mut self, del_names: &[&str]) {
        if del_names.is_empty() || del_names.len() == self.nsamples() {
            panic!("Invalid number of samples to remove")
        }

        // Find position of names in the array rows
        let mut del_name_set = HashSet::new();
        for name in del_names {
            del_name_set.insert(name.to_string());
        }
        let mut idx_list = Vec::new();
        let mut new_names = Vec::new();
        for (idx, name) in self.names.iter_mut().enumerate() {
            if del_name_set.contains(name) {
                idx_list.push(idx);
                del_name_set.remove(name);
            } else {
                new_names.push(mem::take(name));
            }
        }

        if !del_name_set.is_empty() {
            panic!("Could not find sample(s): {:?}", del_name_set);
        }

        let mut idx_it = idx_list.iter();
        let mut next_idx = idx_it.next();
        let mut filtered_variants = Array2::zeros((self.ksize(), 0));
        for (sample_idx, sample_variants) in self.variants.t().outer_iter().enumerate() {
            if let Some(next_idx_val) = next_idx {
                if *next_idx_val == sample_idx {
                    next_idx = idx_it.next();
                    continue;
                }
            }
            filtered_variants.push_column(sample_variants).unwrap();
        }
        self.variants = filtered_variants;
        self.names = new_names;
        self.update_counts(false);
    }

    /// Filters variants (middle bases) by frequency.
    ///
    /// Passing variants (rows) are added to a new array, which overwrites the old one.
    /// Returns the number of removed sites
    ///
    /// # Arguments
    ///
    /// - `min_count` -- minimum number of samples split k-mer found in.
    /// - `filter` -- either none, remove constant, remove ambiguous, or both. See [`FilterType`]
    /// - `mask_ambig` -- replace any non-ACGTUN- with N
    /// - `ignore_const_gaps` -- filter any sites where the only variants are gaps
    /// - `update_kmers` -- update counts and split k-mers after removing variants.
    ///
    /// The default for `update_kmers` should be `true`, but it can be `false`
    /// if [`write_fasta()`] will be called immediately afterwards.
    #[cfg(not(target_arch = "wasm32"))]
    pub fn filter(
        &mut self,
        min_count: usize,
        filter_ambig_as_missing: bool,
        filter: &FilterType,
        mask_ambig: bool,
        ignore_const_gaps: bool,
        update_kmers: bool,
    ) -> i32 {
        if ignore_const_gaps && matches!(filter, FilterType::NoAmbig | FilterType::NoFilter) {
            log::warn!("--no-gap-only-sites can only be applied when filtering constant bases")
        }

        let total = self.names.len();
        let mut filtered_variants = Array2::zeros((0, total));
        let mut filtered_counts = Vec::new();
        let mut filtered_kmers = Vec::new();
        let mut removed = 0;

        if filter_ambig_as_missing {
            self.update_counts(true);
        }

        for count_it in self
            .variant_count
            .iter()
            .zip(self.variants.axis_iter(Axis(0)))
            .zip(self.split_kmers.iter())
        {
            let ((count, row), kmer) = count_it;
            if *count >= min_count {
                let keep_var = match *filter {
                    FilterType::NoFilter => true,
                    FilterType::NoConst => {
                        let mut var_types = HashSet::new();
                        for var in row {
                            if !ignore_const_gaps || *var != b'-' {
                                var_types.insert(*var);
                                if var_types.len() > 1 {
                                    break;
                                }
                            }
                        }
                        var_types.len() > 1
                    }
                    FilterType::NoAmbig => {
                        let mut keep = true;
                        for var in row {
                            if is_ambiguous(*var) {
                                keep = false;
                                break;
                            }
                        }
                        keep
                    }
                    FilterType::NoAmbigOrConst => {
                        let mut var_types = HashSet::new();
                        for var in row {
                            var_types.insert(*var);
                        }
                        let mut count = 0;
                        for base in var_types {
                            let lower_base = base | 0x20;
                            count += match lower_base {
                                b'a' | b'c' | b'g' | b't' | b'u' => 1,
                                b'-' => {
                                    if ignore_const_gaps {
                                        0
                                    } else {
                                        1
                                    }
                                }
                                _ => 0,
                            }
                        }
                        count > 1
                    }
                };
                if keep_var {
                    filtered_variants.push_row(row).unwrap();
                    filtered_counts.push(*count);
                    if update_kmers {
                        filtered_kmers.push(*kmer);
                    }
                } else {
                    removed += 1;
                }
            } else {
                removed += 1;
            }
        }
        self.variants = filtered_variants;
        self.variant_count = filtered_counts;
        if update_kmers {
            self.split_kmers = filtered_kmers;
        }
        log::info!("Filtering removed {removed} split k-mers");

        // Replace any ambiguous variants with Ns, if requested
        if mask_ambig {
            let mut masked = 0;
            self.variants.mapv_inplace(|v| {
                if is_ambiguous(v) {
                    masked += 1;
                    b'N'
                } else {
                    v
                }
            });
            log::info!("Masked {masked} ambiguous bases (non-A/C/G/T/U/N/-) with 'N'");
        }

        removed
    }

    /// Calculates pairwise distances between samples in the array.
    ///
    /// Returns a [`Vec`] of [`Vec`]. This is the upper triangle of the
    /// distance matrix. Each entry is a tuple: first entry is the number of
    /// SNPs different, second entry is the proportion of mismatching k-mers.
    ///
    /// Used with `ska distance`
    ///
    /// # Arguments
    ///
    /// - `constant` – the number of prefiltered constant bases, used to adjust
    ///   the denominator of mismatch proportion
    pub fn distance(&self, constant: f64, filt_ambig: bool) -> Vec<Vec<VariantDist>> {
        let mut distances: Vec<Vec<VariantDist>> = Vec::new();
        self.variants
            .axis_iter(Axis(1))
            .into_par_iter()
            .progress_count(self.variants.ncols() as u64)
            .enumerate()
            .map(|(i, row)| {
                let mut partial_dists: Vec<VariantDist> =
                    Vec::with_capacity(self.variants.ncols() - (i + 1));
                for j in (i + 1)..self.variants.ncols() {
                    partial_dists.push(Self::variant_dist(
                        &row,
                        &self.variants.index_axis(Axis(1), j),
                        constant,
                        filt_ambig,
                    ));
                }
                partial_dists
            })
            .collect_into_vec(&mut distances);
        distances
    }

    /// Removes (weeds) a given set of split k-mers from the array.
    ///
    /// Split k-mers to be removed must be from a single FASTA file (which
    /// may be a multi-FASTA) generated with [`RefSka::new()`]
    ///
    /// Used with `ska weed`.
    ///
    /// # Arguments
    ///
    /// - `weed_ref` -- a processed reference with split k-mers to remove.
    /// - `reverse` -- only remove k-mers not in the input file.
    ///
    pub fn weed(&mut self, weed_ref: &RefSka<IntT>, reverse: bool) {
        let weed_kmers: HashSet<IntT> = HashSet::from_iter(weed_ref.kmer_iter());

        let mut removed = 0;
        let mut new_sk = Vec::new();
        let mut new_variants = Array2::zeros((0, self.nsamples()));
        let mut new_counts = Vec::new();
        for kmer_it in self
            .split_kmers
            .iter()
            .zip(self.variants.outer_iter())
            .zip(self.variant_count.iter())
        {
            let ((kmer, var_row), count) = kmer_it;
            let kmer_found = weed_kmers.contains(kmer);
            if (!reverse && !kmer_found) || (reverse && kmer_found) {
                new_sk.push(*kmer);
                new_variants.push_row(var_row).unwrap();
                new_counts.push(*count);
            } else {
                removed += 1;
            }
        }
        self.split_kmers = new_sk;
        self.variants = new_variants;
        self.variant_count = new_counts;
        if !reverse {
            log::info!("Removed {} of {} weed k-mers", removed, weed_ref.ksize());
        } else {
            log::info!(
                "Kept {} k-mers using {} reverse weed k-mers",
                self.split_kmers.len(),
                weed_ref.ksize()
            );
        }
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
    #[cfg(not(target_arch = "wasm32"))]
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

    /// Write the middle bases as an alignment (FASTA), but not to a file. For wasm.
    #[cfg(target_arch = "wasm32")]
    pub fn write_fasta(&self) -> String {
        // Do the transpose, but recopy with correct strides
        let var_t = self.variants.t();
        let mut var_t_owned = Array2::zeros(var_t.raw_dim());
        var_t_owned.assign(&var_t);
        let mut outputaln = String::new();
        logw(
            &format!(
                "Starting alignment properly, size vars: {:?}, size names: {:?}",
                var_t.len(),
                self.names.len()
            ),
            None,
        );
        self.names
            .iter()
            .zip(var_t_owned.outer_iter())
            .for_each(|it| {
                let (name, seq_u8) = it;

                outputaln.push_str(&("> ".to_owned() + &(name.to_owned() + "\n")));

                outputaln.push_str(
                    &(str::from_utf8(seq_u8.as_slice().expect("Array conversion error"))
                        .unwrap()
                        .to_owned()
                        + "\n"),
                );
            });
        outputaln
    }

    /// Count of number of k-mers found in each sample
    pub fn n_sample_kmers(&self) -> Vec<i32> {
        self.variants
            .map(|v| if *v != b'-' { 1 } else { 0 })
            .sum_axis(Axis(0))
            .to_vec()
    }

    /// K-mer length used when builiding
    pub fn kmer_len(&self) -> usize {
        self.k
    }

    /// Whether reverse complement was used when building
    pub fn rc(&self) -> bool {
        self.rc
    }

    /// Total number of split k-mers
    pub fn ksize(&self) -> usize {
        self.split_kmers.len()
    }

    /// Number of samples
    pub fn nsamples(&self) -> usize {
        self.variants.ncols()
    }

    /// Sample names
    pub fn names(&self) -> &Vec<String> {
        &self.names
    }

    // Distance between two samples (columns of the array)
    fn variant_dist(
        sample1: &ArrayView<u8, Dim<[usize; 1]>>,
        sample2: &ArrayView<u8, Dim<[usize; 1]>>,
        constant: f64,
        filt_ambig: bool,
    ) -> VariantDist {
        //  ACGT vs different ACGT -> +1
        //  Ambig bases are converted to prob vectors and multiplied
        //  '-' vs anything counts as a mismatch
        let mut distance = 0.0;
        let mut mismatches = 0.0;
        let mut matches = constant; // already counted and filtered some matches across whole alignment
        for (var1, var2) in sample1.iter().zip(sample2) {
            if *var1 == b'-' || *var2 == b'-' {
                if !(*var1 == b'-' && *var2 == b'-') {
                    mismatches += 1.0;
                }
            } else if filt_ambig {
                if !is_ambiguous(*var1) && !is_ambiguous(*var2) {
                    matches += 1.0;
                    if *var1 != *var2 {
                        distance += 1.0;
                    }
                }
            } else {
                let var1_p = base_to_prob(*var1);
                let var2_p = base_to_prob(*var2);
                let overlap: f64 = var1_p.iter().zip(var2_p).map(|(p1, p2)| *p1 * p2).sum();
                if overlap > 0.0 {
                    matches += 1.0;
                }
                distance += 1.0 - overlap;
            }
        }
        let mismatch_prop = if (matches + mismatches) == 0.0 {
            0.0
        } else {
            mismatches / (matches + mismatches)
        };
        VariantDist {
            distance,
            mismatch_prop,
            match_count: matches as usize,
            mismatch_count: mismatches as usize,
        }
    }

    /// Iterator over split k-mers and middle bases
    pub fn iter(&self) -> KmerIter<'_, IntT> {
        KmerIter {
            kmers: &self.split_kmers,
            vars: &self.variants,
            index: 0,
        }
    }
}

/// Writes basic information.
///
/// k-mer length, reverse complement, number of split k-mers, number of samples
///
/// Used with `ska nk`
impl<IntT> fmt::Display for MergeSkaArray<IntT>
where
    IntT: for<'a> UInt<'a>,
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "ska_version={}\nk={}\nk_bits={}\nrc={}\nk-mers={}\nsamples={}\n",
            self.ska_version,
            self.kmer_len(),
            self.k_bits,
            self.rc(),
            self.ksize(),
            self.nsamples()
        )?;
        writeln!(f, "sample_names={:?}", self.names)?;
        writeln!(f, "sample_kmers={:?}", self.n_sample_kmers())
    }
}

/// Writes all decoded split k-mers and middle bases.
///
/// Used with `ska nk --full-info`
impl<IntT> fmt::Debug for MergeSkaArray<IntT>
where
    IntT: for<'a> UInt<'a>,
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let (lower_mask, upper_mask) = IntT::generate_masks(self.k);
        self.split_kmers
            .iter()
            .zip(self.variants.outer_iter())
            .try_for_each(|it| {
                let (split_kmer, vars_u8) = it;
                let mut seq_string = String::with_capacity(self.nsamples() * 2);
                for middle_base in vars_u8 {
                    let base = if *middle_base == 0 {
                        '-'
                    } else {
                        *middle_base as char
                    };
                    seq_string.push(base);
                    seq_string.push(',');
                }
                seq_string.pop();
                let (upper, lower) = decode_kmer(self.k, *split_kmer, upper_mask, lower_mask);
                writeln!(f, "{upper}\t{lower}\t{seq_string}")
            })
    }
}

/// Iterator type over split k-mers and middle bases
///
/// Each return is a tuple of the encoded split-kmer and a vector
/// of the encoded middle bases
pub struct KmerIter<'a, IntT> {
    kmers: &'a Vec<IntT>,
    vars: &'a Array2<u8>,
    index: usize,
}

impl<IntT> Iterator for KmerIter<'_, IntT>
where
    IntT: for<'b> UInt<'b>,
{
    // Note this returns a Vec of the middle bases rather than an array
    // because this is more likely to be useful in user code
    type Item = (IntT, Vec<u8>);

    fn next(&mut self) -> Option<Self::Item> {
        if self.index < self.kmers.len() {
            let row = Some((
                self.kmers[self.index],
                self.vars.index_axis(Axis(0), self.index).to_vec(),
            ));
            self.index += 1;
            row
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*; // Import functions and types from the parent module
    use ndarray::array;

    fn setup_struct() -> MergeSkaArray<u64> {
        // Populate with some initial data.
        MergeSkaArray::<u64> {
            k: 31,
            rc: true,
            names: vec![
                "Sample1".to_string(),
                "Sample2".to_string(),
                "Sample3".to_string(),
            ],
            split_kmers: vec![1, 2, 3],
            variants: array![[b'A', b'G', b'Y'], [b'T', b'-', b'Y'], [b'N', b'Y', b'Y']],
            variant_count: vec![3, 2, 3],
            ska_version: "NA".to_string(),
            k_bits: 64,
        }
    }

    #[test]
    fn test_update_counts() {
        let mut merge_ska_array = setup_struct();

        // Test with filter_ambig_as_missing = true
        merge_ska_array.update_counts(true);

        // Check if variant counts are updated correctly
        assert_eq!(merge_ska_array.variant_count, vec![2, 1]); // Expected counts

        // Check if variants are updated correctly
        // Assuming `variants` should now contain only the non-empty rows
        assert_eq!(
            merge_ska_array.variants,
            array![[b'A', b'G', b'Y'], [b'T', b'-', b'Y']]
        );

        // Check if split_kmers are updated correctly
        assert_eq!(merge_ska_array.split_kmers, vec![1, 2]); // Expected kmers
    }

    #[test]
    fn test_kmer_iterator() {
        let ska_array = setup_struct();
        let mut iter = ska_array.iter();

        // First iteration
        let (kmer, vars) = iter.next().unwrap();
        assert_eq!(kmer, 1);
        assert_eq!(vars, vec![b'A', b'G', b'Y']);

        // Second iteration
        let (kmer, vars) = iter.next().unwrap();
        assert_eq!(kmer, 2);
        assert_eq!(vars, vec![b'T', b'-', b'Y']);

        // Third iteration
        let (kmer, vars) = iter.next().unwrap();
        assert_eq!(kmer, 3);
        assert_eq!(vars, vec![b'N', b'Y', b'Y']);

        // No more items
        assert!(iter.next().is_none());
    }

    #[test]
    fn test_delete_samples_normal() {
        let mut ska_array = setup_struct();

        ska_array.delete_samples(&["Sample1", "Sample2"]);

        // Check that the samples were deleted
        assert_eq!(ska_array.names, vec!["Sample3"]);
        assert_eq!(ska_array.variants, array![[b'Y'], [b'Y'], [b'Y']]);
    }

    #[test]
    #[should_panic(expected = "Invalid number of samples to remove")]
    fn test_delete_samples_empty_or_all() {
        let mut ska_array = setup_struct();

        // This should panic
        ska_array.delete_samples(&[]);
    }

    #[test]
    #[should_panic(expected = "Could not find sample(s): ")]
    fn test_delete_samples_non_existent() {
        let mut ska_array = setup_struct();

        // This should panic because "Sample4" does not exist
        ska_array.delete_samples(&["Sample4"]);
    }
}
