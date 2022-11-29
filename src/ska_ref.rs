// Class for split-kmers from one input file
// Used for ska map

use std::fmt;
use std::str;

extern crate needletail;
use needletail::parse_fastx_file;
use ndarray::{ArrayView, Array2};

use crate::merge_ska_dict::MergeSkaDict;
use crate::ska_dict::split_kmer::SplitKmer;
use crate::ska_dict::bit_encoding::RC_IUPAC;

pub struct RefKmer {
    pub kmer: u64,
    pub base: u8,
    pub pos: usize,
    pub chrom: usize,
    pub rc: bool,
}

pub struct RefSka {
    k: usize,
    chrom_names: Vec<String>,
    split_kmer_pos: Vec<RefKmer>, // (kmer, base, pos, rc)
    mapped_variants: Array2<u8>
}

impl RefSka {
    pub fn is_mapped(&self) -> bool {
        self.mapped_variants.nrows() > 0
    }

    pub fn ksize(&self) -> usize {
        self.split_kmer_pos.len()
    }

    // TODO: set some defaults?
    pub fn new(k: usize, filename: &str, rc: bool) -> Self {
        if k < 5 || k > 31 || k % 2 == 0 {
            panic!("Invalid k-mer length");
        }
        let mapped_variants = Array2::zeros((0, 0));
        let mut reader =
            parse_fastx_file(&filename).expect(&format!("Invalid path/file: {}", filename));
        let split_kmer_pos = Vec::new();
        let chrom_names = Vec::new();
        let mut ref_ska = Self {
            k,
            chrom_names,
            split_kmer_pos,
            mapped_variants
        };

        let mut chrom = 0;
        while let Some(record) = reader.next() {
            let seqrec = record.expect("Invalid FASTA record");
            chrom_names.push(str::from_utf8(seqrec.id()).unwrap().to_owned());
            split_kmer_pos.reserve(seqrec.num_bases());
            let kmer_opt = SplitKmer::new(seqrec.seq(), seqrec.num_bases(), k, rc);
            if kmer_opt.is_some() {
                let mut kmer_it = kmer_opt.unwrap();
                let (kmer, base, rc) = kmer_it.get_curr_kmer();
                let mut pos = kmer_it.get_pos();
                split_kmer_pos.push(RefKmer{kmer, base, pos, chrom, rc});
                while let Some((kmer, base, rc)) = kmer_it.get_next_kmer() {
                    pos = kmer_it.get_pos();
                    split_kmer_pos.push(RefKmer{kmer, base, pos, chrom, rc});
                }
            }
            chrom += 1;
        }
        if ref_ska.ksize() == 0 {
            panic!("{} has no valid sequence", filename);
        }
        return ref_ska;
    }

    pub fn map(&mut self, ska_dict: &MergeSkaDict) {
        self.mapped_variants = Array2::zeros((0, ska_dict.nsamples()));
        for ref_k in self.split_kmer_pos {
            if ska_dict.kmer_dict().contains_key(&ref_k.kmer) {
                let seq_char: Vec::<u8> = ska_dict.kmer_dict()[&ref_k.kmer]
                    .iter()
                    .map(|x| {
                        match ref_k.rc {
                            true => RC_IUPAC[*x as usize],
                            false => *x
                        }
                    })
                    .collect();
                self.mapped_variants.push_row(ArrayView::from(&seq_char)).unwrap();
            }
        }
    }

    pub fn write_vcf(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if !self.is_mapped() {
            panic!("Tried to write VCF before variants mapped");
        }
        write!(f, "placeholder")
    }
}

impl fmt::Display for RefSka {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if !self.is_mapped() {
            panic!("Tried to write VCF before variants mapped");
        }
        // TODO: write out by editing variants into each ref row
        write!(f, "placeholder")
    }
}
