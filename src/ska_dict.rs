// Class for split-kmers from one input file

extern crate needletail;
use needletail::parse_fastx_file;
use hashbrown::HashMap;

pub mod split_kmer;
use crate::ska_dict::split_kmer::SplitKmer;

pub mod bit_encoding;
use crate::ska_dict::bit_encoding::{encode_base, IUPAC};

pub struct SkaDict {
    sample_idx: usize,
    name: String,
    split_kmers: HashMap<u64, u8>
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

    pub fn idx(&self) -> usize {
        self.sample_idx
    }

    pub fn kmers(&self) -> &HashMap<u64, u8> {
        &self.split_kmers
    }

    pub fn name(&self) -> &String {
        &self.name
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
