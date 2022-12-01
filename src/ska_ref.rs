// Class for split-kmers from one input file
// Used for ska map

use std::fmt;
use std::str;

extern crate needletail;
use needletail::parse_fastx_file;
use ndarray::{ArrayView, Array2, s};

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
    // Index
    k: usize,
    split_kmer_pos: Vec<RefKmer>,

    // Sequence
    chrom_names: Vec<String>,
    seq: Vec<Vec<u8>>,
    total_size: usize,

    // Mapping
    mapped_pos: Vec<(usize, usize)>, // (chrom, pos)
    mapped_variants: Array2<u8>,
    mapped_names: Vec<String>
}

// TODO want a mapped ska
pub struct MapSka {

}

impl RefSka {
    fn is_mapped(&self) -> bool {
        self.mapped_variants.nrows() > 0
    }

    pub fn ksize(&self) -> usize {
        self.split_kmer_pos.len()
    }

    pub fn new(k: usize, filename: &str, rc: bool) -> Self {
        if k < 5 || k > 31 || k % 2 == 0 {
            panic!("Invalid k-mer length");
        }

        let split_kmer_pos = Vec::new();
        let seq = Vec::new();
        let chrom_names = Vec::new();
        let mut total_size = 0;

        let mut reader =
            parse_fastx_file(&filename).expect(&format!("Invalid path/file: {}", filename));
        let mut chrom = 0;
        while let Some(record) = reader.next() {
            let seqrec = record.expect("Invalid FASTA record");
            chrom_names.push(str::from_utf8(seqrec.id()).unwrap().to_owned());
            split_kmer_pos.reserve(seqrec.num_bases());
            total_size += seqrec.num_bases();

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
            seq.push(seqrec.seq().to_vec());
        }
        if split_kmer_pos.len() == 0 {
            panic!("{} has no valid sequence", filename);
        }

        let mapped_variants = Array2::zeros((0, 0));
        let mapped_pos = Vec::new();
        let mapped_names = Vec::new();
        Self {
            k,
            seq,
            total_size,
            chrom_names,
            split_kmer_pos,
            mapped_pos,
            mapped_variants,
            mapped_names
        }
    }

    pub fn map(&mut self, ska_dict: &MergeSkaDict) {
        self.mapped_names = *ska_dict.names();
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
                self.mapped_pos.push((ref_k.chrom, ref_k.pos));
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
        if self.chrom_names.len() > 1 {
            eprintln!("WARNING: Reference contained multiple contigs, in the output they will be concatenated");
        }
        for sample_name in self.mapped_names {
            let sample_vars = self.mapped_variants.slice(s![.., sample_idx]);
            let seq: Vec<u8> = Vec::new();
            seq.reserve(self.total_size);

            // TODO fix for multi-chrom case
            let (mut next_pos, map_pos) = (0, 0);
            self.mapped_pos.iter()
                .zip(sample_vars.iter())
                .for_each(|it| {
                    let ((map_chrom, map_pos), base) = it;
                    if *map_pos > next_pos {
                        // Copy in ref seq if no k-mers mapped over a region
                        seq.extend_from_slice(&self.seq[0][next_pos..*map_pos]);
                    }
                    next_pos = *map_pos + 1;
                    seq.push(*base);
                });
            let seq_char: String = seq.iter().map(|x| *x as char).collect();
            write!(f, ">{}\n{}\n", sample_name, seq_char);
        }
    }
}
