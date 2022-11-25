// Class for split-kmers from one input file

extern crate needletail;
use hashbrown::HashMap;
use needletail::parse_fastx_file;

pub mod split_kmer;
use crate::ska_dict::split_kmer::SplitKmer;

pub mod bit_encoding;
use crate::ska_dict::bit_encoding::{decode_base, IUPAC};

pub struct SkaDict {
    k: usize,
    sample_idx: usize,
    name: String,
    split_kmers: HashMap<u64, u8>,
    positions: Option<HashMap<u64, usize>>,
}

impl SkaDict {
    fn add_to_dict(&mut self, kmer: u64, base: u8, pos: usize) {
        self.split_kmers
            .entry(kmer)
            .and_modify(|b| *b = IUPAC[base as usize * 256 + *b as usize])
            .or_insert(decode_base(base));
        if self.positions.is_some() {
            self.positions.unwrap().entry(kmer).or_insert(pos);
        }
    }

    pub fn kmer_len(&self) -> usize {
        self.k
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

    // TODO: set some defaults?
    pub fn new(k: usize, sample_idx: usize, filename: &str, name: &str, rc: bool, keep_pos: bool) -> Self {
        if k < 5 || k > 31 || k % 2 == 0 {
            panic!("Invalid k-mer length");
        }
        let mut reader =
            parse_fastx_file(&filename).expect(&format!("Invalid path/file: {}", filename));
        let name = name.to_string();
        let split_kmers: HashMap<u64, u8> = HashMap::default();
        let mut positions = None;
        if keep_pos {
            positions = Some(HashMap::default());
        }
        let mut sk_dict = Self {
            k,
            sample_idx,
            name,
            split_kmers,
            positions,
        };
        while let Some(record) = reader.next() {
            let seqrec = record.expect("Invalid FASTA record");
            let kmer_opt = SplitKmer::new(seqrec.seq(), seqrec.num_bases(), k, rc);
            if kmer_opt.is_some() {
                let mut kmer_it = kmer_opt.unwrap();
                let (kmer, base) = kmer_it.get_curr_kmer();
                sk_dict.add_to_dict(kmer, base, kmer_opt.get_pos());
                while let Some((kmer, base)) = kmer_it.get_next_kmer() {
                    sk_dict.add_to_dict(kmer, base, kmer_opt.get_pos());
                }
            }
        }
        if sk_dict.ksize() == 0 {
            panic!("{} has no valid sequence", filename);
        }
        return sk_dict;
    }
}
