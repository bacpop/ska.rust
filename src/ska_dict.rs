
use std::fmt;
use std::mem;
use core::mem::swap;

use std::error::Error;
use std::fs::File;
use std::io::{BufWriter, BufReader};

extern crate needletail;
use needletail::parse_fastx_file;
use hashbrown::HashMap;
use ndarray::{Array2, ArrayView, Axis};
use serde::{Serialize, Deserialize};

pub mod split_kmer;
use crate::ska_dict::split_kmer::SplitKmer;

pub mod bit_encoding;
use crate::ska_dict::bit_encoding::{encode_base, decode_kmer, IUPAC};

pub struct SkaDict {
    sample_idx: usize,
    name: String,
    split_kmers: HashMap<u64, u8>
}

pub struct MergeSkaDict {
    n_samples: usize,
    names: Vec<String>,
    split_kmers: HashMap<u64, Vec<u8>>
}

#[derive(Serialize, Deserialize, Debug)]
pub struct MergeSkaArray {
    names: Vec<String>,
    split_kmers: Vec<u64>,
    variants: Array2<u8>,
    variant_count: Vec<usize>
}


impl SkaDict {
    fn add_to_dict(&mut self, kmer: u64, base: u8) {
        self.split_kmers.entry(kmer)
            .and_modify(|b| {*b = IUPAC[encode_base(*b) as usize * 256 + base as usize]})
            .or_insert(base);
    }

    pub fn ksize(&self) -> usize {
        self.split_kmers.len()
    }

    pub fn new(sample_idx: usize, filename: &str, name: &str, rc: bool) -> Self {
        let mut reader = parse_fastx_file(&filename).expect(&format!("Invalid path/file: {}", filename));
        let name = name.to_string();
        let split_kmers: HashMap<u64, u8> = HashMap::default();
        let mut sk_dict = Self{sample_idx, name, split_kmers};
        while let Some(record) = reader.next() {
            let seqrec = record.expect("Invalid FASTA record");
            let kmer_opt = SplitKmer::new(seqrec.seq(), seqrec.num_bases(), rc);
            if kmer_opt.is_some() {
                let mut kmer_it = kmer_opt.unwrap();
                let (kmer, base) = kmer_it.get_curr_kmer();
                sk_dict.add_to_dict(kmer, base);
                while let Some((kmer, base)) = kmer_it.get_next_kmer() {
                    sk_dict.add_to_dict(kmer, base);
                }
            }
        }
        if sk_dict.ksize() == 0 {
            panic!("{} has no valid sequence", filename);
        }
        return sk_dict;
    }
}

impl MergeSkaDict {
    pub fn ksize(&self) -> usize {
        self.split_kmers.len()
    }

    pub fn nsamples(&self) -> usize {
        self.n_samples
    }

    pub fn new(n_samples: usize) -> Self {
        let names = vec!["".to_string(); n_samples];
        let split_kmers = HashMap::default();
        return Self {n_samples, names, split_kmers};
    }

    pub fn append(&mut self, other: &SkaDict) {
        self.names[other.sample_idx] = other.name.clone();
        if self.ksize() == 0 {
            for (kmer, base) in &other.split_kmers {
                let mut base_vec: Vec<u8> = vec![b'-'; self.n_samples];
                base_vec[other.sample_idx] = *base;
                self.split_kmers.insert(*kmer, base_vec);
            }
        } else {
            for (kmer, base) in &other.split_kmers {
                self.split_kmers.entry(*kmer)
                    .and_modify(|b| {
                        b[other.sample_idx] = *base;
                    })
                    .or_insert_with(|| {
                        let mut new_base_vec: Vec<u8> = vec![b'-'; self.n_samples];
                        new_base_vec[other.sample_idx] = *base;
                        new_base_vec });
            }
        }
    }

    pub fn merge<'a>(&'a mut self, other: &'a mut MergeSkaDict) {
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
                    self.split_kmers.entry(*kmer)
                        .and_modify(|self_vec| {
                            for base_it in other_vec.iter().zip(self_vec.iter_mut()) {
                                let (other_base, self_base) = base_it;
                                if *self_base == b'N' {
                                    *self_base = *other_base;
                                }
                            }
                        })
                        .or_insert(mem::take(other_vec));
                }
            }
        }
    }
}

// TODO: take a formatter which determines whether k-mers are written out
impl fmt::Display for MergeSkaDict {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}\tkmers\n{}\tsamples\n",
                   self.ksize(), self.nsamples())?;
        writeln!(f, "{:?}", self.names)?;
        self.split_kmers.iter().try_for_each(|it| {
            let (split_kmer, vars_u8) = it;
            let mut seq_string = String::with_capacity(self.nsamples());
            for middle_base in vars_u8 {
                seq_string.push(*middle_base as char);
                seq_string.push(',');
            }
            seq_string.pop();
            let (upper, lower) = decode_kmer(*split_kmer);
            write!(f, "{}\t{}\t{}\n", upper, lower, seq_string)
        })
    }
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
        let names = dynamic.names.clone();
        let mut variants = Array2::zeros((0, dynamic.nsamples()));
        let mut split_kmers: Vec<u64> = Vec::new();
        split_kmers.reserve(dynamic.ksize());
        let mut variant_count: Vec<usize> = Vec::new();
        // TODO parallel it -> use slice_mut
        for (kmer, bases) in &dynamic.split_kmers {
            split_kmers.push(*kmer);
            variant_count.push(bases.iter().filter(|b| **b != b'-').count());
            variants.push_row(ArrayView::from(bases)).unwrap();

        }
        Self {names, split_kmers, variants, variant_count}
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
        let names = self.names.clone();
        let mut split_kmers: HashMap<u64, Vec<u8>> = HashMap::new();
        split_kmers.reserve(self.variants.nrows());
        //TODO parallelise
        for row_it in self.variants.outer_iter().zip(self.split_kmers.iter()) {
            let (row_vec, kmer) = row_it;
            split_kmers.insert(*kmer, row_vec.to_vec());
        }
        MergeSkaDict{n_samples, names, split_kmers}
    }

    pub fn delete_samples(&mut self, idx: &mut Vec<usize>) {
        idx.sort();
        let mut idx_it = idx.iter();
        let mut next_idx = idx_it.next();

        let new_size = self.names.len() - idx.len();
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

    pub fn map(&self, ref_fasta: &str) {

    }

    pub fn filter(&mut self, min_count: usize, const_sites: bool) {
        let total = self.names.len();
        let mut filtered_variants = Array2::zeros((0, total));
        let mut filtered_counts = Vec::new();
        for count_it in self.variant_count.iter().zip(self.variants.axis_iter(Axis(0))) {
            let (count, row) = count_it;
            if *count >= min_count {
                if !const_sites {
                    let first_var = row[0];
                    for var in row {
                        if *var != first_var {
                            break;
                        }
                    }
                }
                filtered_variants.push_row(row).unwrap();
                filtered_counts.push(*count);
            }
        }
        self.variants = filtered_variants;
        self.variant_count = filtered_counts;
    }
}

impl fmt::Display for MergeSkaArray {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        self.names.iter()
            .zip(self.variants.t().outer_iter())
            .try_for_each(|it| {
                let (name, seq_u8) = it;
                let seq_char: String = seq_u8.iter().map(|x| *x as char).collect();
                write!(f, ">{}\n{}\n", name, seq_char)
            })
    }
}
