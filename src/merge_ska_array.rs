// Class for split-kmers from multiple SkaDict
// In array representation to support:
//  filter
//  save/load
//  mapping
//  sample deletion
// Can be converted to/from MergeSkaDict
// Print will print out alignment

use std::error::Error;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};

use hashbrown::HashMap;
use ndarray::{Array2, ArrayView, Axis};
use needletail::parser::write_fasta;
use serde::{Deserialize, Serialize};

use crate::merge_ska_dict::MergeSkaDict;

#[derive(Serialize, Deserialize, Debug)]
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
        let mut new_variants = Array2::zeros((0, self.names.len()));
        for var_row in self.variants.outer_iter() {
            let count = var_row.iter().filter(|b| **b != b'-').count();
            if count > 0 {
                new_counts.push(count);
                new_variants.push_row(var_row).unwrap();
            }
        }
        self.variants = new_variants;
        self.variant_count = new_counts;
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
        // TODO parallel it -> use slice_mut
        for (kmer, bases) in dynamic.kmer_dict() {
            split_kmers.push(*kmer);
            variant_count.push(bases.iter().filter(|b| **b != b'-').count());
            variants.push_row(ArrayView::from(bases)).unwrap();
        }
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
        let mut serial_file = BufWriter::new(File::create(filename)?);
        ciborium::ser::into_writer(self, &mut serial_file)?;
        Ok(())
    }

    pub fn load(filename: &str) -> Result<Self, Box<dyn Error>> {
        let ska_file = BufReader::new(File::open(filename)?);
        let ska_obj: Self = ciborium::de::from_reader(ska_file)?;
        Ok(ska_obj)
    }

    pub fn to_dict(&self) -> MergeSkaDict {
        let n_samples = self.names.len();
        let mut names = self.names.clone();
        let mut split_kmers: HashMap<u64, Vec<u8>> = HashMap::new();
        split_kmers.reserve(self.variants.nrows());
        //TODO parallelise
        for row_it in self.variants.outer_iter().zip(self.split_kmers.iter()) {
            let (row_vec, kmer) = row_it;
            split_kmers.insert(*kmer, row_vec.to_vec());
        }
        let mut dict = MergeSkaDict::new(self.k, n_samples, self.rc);
        dict.build_from_array(&mut names, &mut split_kmers);
        return dict;
    }

    pub fn delete_samples(&mut self, del_names: &Vec<String>) {
        // Find position of names in the array rows
        let mut name_dict: HashMap<String, usize> = HashMap::default();
        for (idx, name_pair) in self.names.iter().enumerate() {
            name_dict.insert(name_pair.to_string(), idx);
        }

        let mut idx_list = Vec::new();
        idx_list.reserve(del_names.len());
        for name in del_names {
            match name_dict.get(name) {
                Some(idx) => idx_list.push(*idx),
                None => panic!("Could not find sample {name}"),
            }
        }
        idx_list.sort();
        let mut idx_it = idx_list.iter();
        let mut next_idx = idx_it.next();

        let new_size = self.names.len() - idx_list.len();
        let mut filtered_variants = Array2::zeros((0, new_size));
        for sample_it in self.variants.t().outer_iter().enumerate() {
            let (s_idx, sample_variants) = sample_it;
            if *next_idx.unwrap_or(&new_size) == s_idx {
                filtered_variants.push_column(sample_variants).unwrap();
            } else {
                next_idx = idx_it.next();
            }
        }
        self.variants = filtered_variants;
        self.update_counts();
    }

    pub fn filter(&mut self, min_count: usize, const_sites: bool) {
        let total = self.names.len();
        let mut filtered_variants = Array2::zeros((0, total));
        let mut filtered_counts = Vec::new();
        for count_it in self
            .variant_count
            .iter()
            .zip(self.variants.axis_iter(Axis(0)))
        {
            let (count, row) = count_it;
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
            }
        }
        self.variants = filtered_variants;
        self.variant_count = filtered_counts;
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
