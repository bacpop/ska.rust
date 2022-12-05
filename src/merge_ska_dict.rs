// Class for split-kmers from multiple SkaDict
// In dictonary representation to support merging/adding
// Print will print out no kmers, no samples, split k-mers

use core::mem::swap;
use core::panic;
use std::fmt;
use std::mem;

use hashbrown::HashMap;

use crate::ska_dict::bit_encoding::{generate_masks, decode_kmer};
use crate::ska_dict::SkaDict;

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
                let mut base_vec: Vec<u8> = vec![b'-'; self.n_samples];
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
                        let mut new_base_vec: Vec<u8> = vec![b'-'; self.n_samples];
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
                            for base_it in other_vec.iter().zip(self_vec.iter_mut()) {
                                let (other_base, self_base) = base_it;
                                if *self_base == b'-' {
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
                seq_string.push(*middle_base as char);
                seq_string.push(',');
            }
            seq_string.pop();
            let (upper, lower) = decode_kmer(self.k, *split_kmer, upper_mask, lower_mask);
            write!(f, "{}\t{}\t{}\n", upper, lower, seq_string)
        })
    }
}
