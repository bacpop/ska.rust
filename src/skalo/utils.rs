//! Data structures used in ska lo
use std::path::PathBuf;
use std::sync::OnceLock;

use crate::ska_dict::bit_encoding::{decode_base, encode_base};

/// Structure to store parameters for skalo
#[derive(Debug)]
pub struct Config {
    /// Input file name
    pub input_file: String,
    /// Output file name
    pub output_name: String,
    /// Maximum missing data graction
    pub max_missing: f32,
    /// Maximum depth of recursive paths
    pub max_depth: usize,
    /// Maximum number of internal indel k-mers
    pub max_indel_kmers: usize,
    /// Number of CPU threads
    pub nb_threads: usize,
    /// Reference genome to map against
    pub reference_genome: Option<PathBuf>,
}
/// Single [`Config`]
pub static CONFIG: OnceLock<Config> = OnceLock::new();

/// Structure to store dataset information
#[derive(Debug, Clone)]
pub struct DataInfo {
    /// k used for graph
    pub k_graph: usize,
    /// Sample names in graph
    pub sample_names: Vec<String>,
}
/// Single [`DataInfo`]
pub static DATA_INFO: OnceLock<DataInfo> = OnceLock::new();

/// Returns the reverse-complement of a DNA string slice
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

/// Structure to store variant groups
#[derive(Clone)]
pub struct VariantInfo {
    /// Reference sequence
    pub sequence: DnaSequence,
    /// SNP postitions
    pub vec_snps: Vec<usize>,
}

impl VariantInfo {
    /// Create a struct holding a reference sequence and SNPs against it
    pub fn new(sequence: DnaSequence, vec_snps: Vec<usize>) -> Self {
        VariantInfo { sequence, vec_snps }
    }
}

/// structure to store DNA sequence in a bit-packed [u8]
#[derive(Clone)]
pub struct DnaSequence {
    /// Sequence as bytes
    pub data: Vec<u8>,
    /// Length of bases
    pub original_length: usize,
}

impl DnaSequence {
    /// Create a new DnaSequence from a DNA string
    pub fn encode(dna: &str) -> Self {
        let data: Vec<u8> = dna.as_bytes().iter().map(|nt| encode_base(*nt)).collect();

        DnaSequence {
            data,
            original_length: dna.len(),
        }
    }

    /// Decode the bit-packed DnaSequence back into a DNA string
    pub fn decode(&self) -> String {
        let out: String = self
            .data
            .iter()
            .map(|nt| decode_base(*nt & 0b11) as char)
            .collect();
        out
    }

    /// Return the original number of nucleotides
    pub fn len(&self) -> usize {
        self.original_length
    }

    /// If the sequence is empty
    pub fn is_empty(&self) -> bool {
        self.original_length == 0
    }

    /// Get the decoded sequence in the provided interval
    pub fn get_range(&self, start: usize, end: usize) -> Vec<u8> {
        let out: Vec<u8> = self.data[start..end]
            .iter()
            .map(|nt| decode_base(*nt & 0b11))
            .collect();
        out
    }
}
