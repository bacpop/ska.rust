
use std::fmt;
extern crate needletail;
use needletail::{parse_fastx_file, Sequence};
use hashbrown::HashMap;

pub mod split_kmer;
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
            .and_modify(|b| {*b[0] = iupac[encode_base(b[0]) * 256 + base]})
            .or_insert(Vec::from([base]));
    }

    pub fn new(&filename: str, &name: str, rc: bool) -> Self {
        let mut reader = parse_fastx_file(&filename).expect(format!("invalid path/file: {}", filename));
        let names = Vec::from([name]);
        let mut sk_dict: HashMap<u64, Vec<u8>> = HashMap::new();
        while let Some(record) = reader.next() {
            let seqrec = record.expect("Invalid FASTA record");
            let kmer_it = SplitKmer(seqrec.seq(), seqrec.num_bases(), rc);
            let (mut kmer, mut base) = kmer_it.get_curr_kmer();
            self.add_to_dict(kmer, base);
            while Some(kmer, base) = kmer_it.get_next_kmer() {
                self.add_to_dict(kmer, base);
            }
        }
        return Self{names, sk_dict};
    }

    pub fn merge(&mut self, &other: SkaDict) {

    }

    pub fn filter(&mut self, min_af: i32) {

    }

    pub fn delete_samples(&mut self, &idx: Vec<u32>) {

    }

    pub fn map(&self, &ref_fasta: str) {

    }

    pub fn serialise(&self) {

    }
}

// This will eventually print the alignment
impl fmt::Display for SkaDict {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?}", self.names)
    }
}