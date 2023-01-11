// Class for split-kmers from multiple SkaDict
// In array representation to support:
//  filter
//  save/load
//  mapping
//  sample deletion
// Can be converted to/from MergeSkaDict
// Print will print out alignment

use core::panic;
use std::error::Error;
use std::fmt;
use std::fs::File;
use std::mem;
use std::io::{BufReader, BufWriter, Write};

use hashbrown::{HashMap, HashSet};
use ndarray::{Array2, ArrayView, Axis};
use needletail::parser::write_fasta;
use serde::{Deserialize, Serialize};

use crate::merge_ska_dict::MergeSkaDict;
use crate::ska_dict::bit_encoding::{decode_kmer, generate_masks};
use crate::ska_ref::RefSka;

#[derive(Serialize, Deserialize)]
pub struct MergeSkaArray {
    k: usize,
    rc: bool,
    names: Vec<String>,
    split_kmers: Vec<u64>,
    variants: Array2<u8>,
    variant_count: Vec<usize>,
}

impl MergeSkaArray {
    fn update_counts(&mut self) {
        let mut new_counts = Vec::new();
        new_counts.reserve(self.variant_count.len());

        let mut new_sk = Vec::new();
        new_sk.reserve(self.split_kmers.len());

        let mut new_variants = Array2::zeros((0, self.names.len()));
        for (var_row, sk) in self.variants.outer_iter().zip(self.split_kmers.iter()) {
            let count = var_row.iter().filter(|b| **b != b'-').count();
            if count > 0 {
                new_counts.push(count);
                new_sk.push(*sk);
                new_variants.push_row(var_row).unwrap();
            }
        }
        self.split_kmers = new_sk;
        self.variants = new_variants;
        self.variant_count = new_counts;
    }

    pub fn kmer_len(&self) -> usize {
        self.k
    }

    pub fn rc(&self) -> bool {
        self.rc
    }

    pub fn ksize(&self) -> usize {
        self.split_kmers.len()
    }

    pub fn nsamples(&self) -> usize {
        self.variants.ncols()
    }

    pub fn new(dynamic: &MergeSkaDict) -> Self {
        let k = dynamic.kmer_len();
        let rc = dynamic.rc();
        let names = dynamic.names().clone();
        let mut variants = Array2::zeros((0, dynamic.nsamples()));
        let mut split_kmers: Vec<u64> = Vec::new();
        split_kmers.reserve(dynamic.ksize());
        let mut variant_count: Vec<usize> = Vec::new();
        for (kmer, bases) in dynamic.kmer_dict() {
            split_kmers.push(*kmer);
            variant_count.push(bases.iter().filter(|b| **b != 0 && **b != b'-').count());
            variants.push_row(ArrayView::from(bases)).unwrap();
        }
        variants.mapv_inplace(|b| u8::max(b, b'-')); // turns zeros to missing
        Self {
            k,
            rc,
            names,
            split_kmers,
            variants,
            variant_count,
        }
    }

    pub fn save(&self, filename: &str) -> Result<(), Box<dyn Error>> {
        let serial_file = BufWriter::new(File::create(filename)?);
        let mut compress_writer = snap::write::FrameEncoder::new(serial_file);
        ciborium::ser::into_writer(self, &mut compress_writer)?;
        Ok(())
    }

    pub fn load(filename: &str) -> Result<Self, Box<dyn Error>> {
        let ska_file = BufReader::new(File::open(filename)?);
        let decompress_reader = snap::read::FrameDecoder::new(ska_file);
        let ska_obj: Self = ciborium::de::from_reader(decompress_reader)?;
        Ok(ska_obj)
    }

    pub fn to_dict(&self) -> MergeSkaDict {
        let n_samples = self.names.len();
        let mut names = self.names.clone();
        let mut split_kmers: HashMap<u64, Vec<u8>> = HashMap::new();
        split_kmers.reserve(self.variants.nrows());
        for row_it in self.variants.outer_iter().zip(self.split_kmers.iter()) {
            let (row_vec, kmer) = row_it;
            split_kmers.insert(*kmer, row_vec.to_vec());
        }
        let mut dict = MergeSkaDict::new(self.k, n_samples, self.rc);
        dict.build_from_array(&mut names, &mut split_kmers);
        return dict;
    }

    pub fn delete_samples(&mut self, del_names: &[&str]) {
        if del_names.len() == 0 || del_names.len() == self.nsamples() {
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

        if del_name_set.len() > 0 {
            panic!("Could not find sample(s): {:?}", del_name_set);
        }
        self.names = new_names;

        let mut idx_it = idx_list.iter();
        let mut next_idx = idx_it.next();
        let new_size = self.names.len() - idx_list.len();
        let mut filtered_variants = Array2::zeros((self.ksize(), 0));
        for (sample_idx, sample_variants) in self.variants.t().outer_iter().enumerate() {
            if *next_idx.unwrap_or(&new_size) == sample_idx {
                next_idx = idx_it.next();
            } else {
                filtered_variants.push_column(sample_variants).unwrap();
            }
        }
        self.variants = filtered_variants;
        self.update_counts();
    }

    pub fn filter(&mut self, min_count: usize, const_sites: bool, update_kmers: bool) {
        let total = self.names.len();
        let mut filtered_variants = Array2::zeros((0, total));
        let mut filtered_counts = Vec::new();
        let mut filtered_kmers = Vec::new();
        let mut removed = 0;
        for count_it in self
            .variant_count
            .iter()
            .zip(self.variants.axis_iter(Axis(0)))
            .zip(self.split_kmers.iter())
        {
            let ((count, row), kmer) = count_it;
            let mut keep_var = false;
            if *count >= min_count {
                if !const_sites {
                    let first_var = row[0];
                    for var in row {
                        if *var != first_var {
                            keep_var = true;
                            break;
                        }
                    }
                } else {
                    keep_var = true;
                }
            }
            if keep_var {
                filtered_variants.push_row(row).unwrap();
                filtered_counts.push(*count);
                if update_kmers {
                    filtered_kmers.push(*kmer);
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
    }

    pub fn weed(&mut self, weed_ref: &RefSka) {
        let weed_kmers: HashSet<u64> = HashSet::from_iter(weed_ref.kmer_iter());

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
            if !weed_kmers.contains(kmer) {
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
        log::info!("Removed {} of {} weed k-mers", removed, weed_ref.ksize());
    }

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
}

impl fmt::Display for MergeSkaArray {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "k={}\nrc={}\n{} k-mers\n{} samples\n",
            self.kmer_len(),
            self.rc(),
            self.ksize(),
            self.nsamples()
        )?;
        writeln!(f, "{:?}", self.names)
    }
}

impl fmt::Debug for MergeSkaArray {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let (lower_mask, upper_mask) = generate_masks(self.k);
        println!("{}", self.variants);
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
                write!(f, "{}\t{}\t{}\n", upper, lower, seq_string)
            })
    }
}
