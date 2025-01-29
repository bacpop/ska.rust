use hashbrown::HashMap;
use std::path::PathBuf;
use std::sync::OnceLock;

use crate::ska_dict::bit_encoding::{decode_base, encode_base};

// structure to hold arguments
#[derive(Debug)]
pub struct Config {
    pub input_file: String,
    pub output_name: String,
    pub max_missing: f32,
    pub max_depth: usize,
    pub max_indel_kmers: usize,
    pub nb_threads: usize,
    pub reference_genome: Option<PathBuf>,
}
pub static CONFIG: OnceLock<Config> = OnceLock::new();

// structure to hold dataset information
#[derive(Debug, Clone)]
pub struct DataInfo {
    pub k_graph: usize,
    pub sample_names: Vec<String>,
}
pub static DATA_INFO: OnceLock<DataInfo> = OnceLock::new();

pub fn rev_compl(seq: &str) -> String {
    let out: String = seq
        .chars()
        .rev()
        .map(|nt| match nt {
            'A' => 'T',
            'C' => 'G',
            'T' => 'A',
            'G' => 'C',
            _ => panic!("Error taking reverse complement of {}", nt),
        })
        .collect();

    out
}

#[derive(Clone)]
pub struct VariantInfo {
    pub sequence: DnaSequence,
    pub vec_snps: Vec<usize>,
}

impl VariantInfo {
    pub fn new(sequence: DnaSequence, vec_snps: Vec<usize>) -> Self {
        VariantInfo { sequence, vec_snps }
    }
}

/// structure to store DNA sequence in a bit-packed [u8]
#[derive(Clone)]
pub struct DnaSequence {
    pub data: Vec<u8>,
    pub original_length: usize,
}

impl DnaSequence {
    // create a new DnaSequence from a DNA string
    pub fn encode(dna: &str) -> Self {
        let data: Vec<u8> = dna.as_bytes().iter().map(|nt| encode_base(*nt)).collect();

        DnaSequence {
            data,
            original_length: dna.len(),
        }
    }

    // decode the bit-packed DnaSequence back into a DNA string
    pub fn decode(&self) -> String {
        let out: String = self
            .data
            .iter()
            .map(|nt| decode_base(*nt & 0b11) as char)
            .collect();
        out
    }

    // return the original number of nucleotides
    pub fn len(&self) -> usize {
        self.original_length
    }

    pub fn get_range(&self, start: usize, end: usize) -> Vec<u8> {
        let out: Vec<u8> = self.data[start..end]
            .iter()
            .map(|nt| decode_base(*nt & 0b11))
            .collect();
        out
    }
}
