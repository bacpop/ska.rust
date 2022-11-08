
use std::fmt;
use core::mem::swap;

extern crate needletail;
use needletail::parse_fastx_file;
use hashbrown::HashMap;

pub mod split_kmer;
use crate::ska_dict::split_kmer::SplitKmer;

pub mod bit_encoding;
use crate::ska_dict::bit_encoding::encode_base;
use crate::ska_dict::bit_encoding::IUPAC;

// Alternative (may be faster? not sure) have BitVec with n_samples allocated
// placing bit in the correct column when reading in
// Then if found can just OR the bitvecs
// If not found just move/copy the bitvec from the new one

// When writing structure, we can just try the compress module on the full split_kmer dict
// BitVecs should probably be four big arrays: maybe just compress these (could reorder columns)

// This will do for everything
// Also thought about having four vec of bvs, but need to make sure
// order is the same as the k-mers, which also need quick lookup (could be done with an index)

#[derive(Default, Clone)]
pub struct SkaDict {
    names: Vec<String>,
    split_kmers: HashMap<u64, Vec<u8>>
}

impl SkaDict {
    fn add_to_dict(&mut self, kmer: u64, base: u8) {
        self.split_kmers.entry(kmer)
            .and_modify(|b| {b[0] = IUPAC[encode_base(b[0]) as usize * 256 + base as usize]})
            .or_insert(Vec::from([base]));
    }

    pub fn ksize(&self) -> usize {
        self.split_kmers.len()
    }

    pub fn nsamples(&self) -> usize {
        self.names.len()
    }

    pub fn new(filename: &str, name: &str, rc: bool) -> Self {
        let mut reader = parse_fastx_file(&filename).expect(&format!("Invalid path/file: {}", filename));
        let names = Vec::from([name.to_string()]);
        let split_kmers: HashMap<u64, Vec<u8>> = HashMap::new();
        let mut sk_dict = Self{names, split_kmers};
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

    pub fn merge<'a>(&'a mut self, other: &'a mut SkaDict) {
        if other.nsamples() > 0 {
            if self.names.is_empty() {
                swap(&mut other.names, &mut self.names);
                swap(&mut other.split_kmers, &mut self.split_kmers);
            } else {
                self.names.append(&mut other.names);
                // NB - if this is changed to a B-Tree need to change this
                // into an ordered merge
                for (kmer, base_vec) in &mut other.split_kmers {
                    self.split_kmers.entry(*kmer)
                        .and_modify(|b| {b.append(base_vec)})
                        .or_insert_with(|| {
                            let mut new_base_vec: Vec<u8> = vec![b'N'; self.names.len()];
                            new_base_vec.append(base_vec);
                            new_base_vec });
                }
                for (kmer, base_vec) in &mut self.split_kmers {
                    if !other.split_kmers.contains_key(kmer) {
                        base_vec.extend(vec![b'N'; other.names.len()]);
                    }
                }
            }
        }
    }

    pub fn filter(&mut self, min_af: i32) {

    }

    pub fn delete_samples(&mut self, idx: &Vec<u32>) {

    }

    pub fn map(&self, ref_fasta: &str) {

    }

    pub fn serialise(&self) {

    }
}

// This will eventually print the alignment
impl fmt::Display for SkaDict {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?}\n", self.names)
    }
}