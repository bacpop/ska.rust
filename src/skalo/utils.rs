use hashbrown::HashMap;
use std::path::PathBuf;
use std::sync::OnceLock;

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
    let reverse_code: HashMap<char, char> =
        [('A', 'T'), ('C', 'G'), ('G', 'C'), ('T', 'A'), ('N', 'N')]
            .iter()
            .cloned()
            .collect();

    let reversed_seq: String = seq.chars().rev().collect();
    let complemented_seq: String = reversed_seq.chars().map(|n| reverse_code[&n]).collect();

    complemented_seq
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
        let mut data = Vec::with_capacity((dna.len() + 3) / 4); // 4 nucleotides per byte
        let mut current_byte = 0u8;
        let mut shift = 6; // start with the highest 2 bits

        for (i, nucleotide) in dna.chars().enumerate() {
            let bits = match nucleotide {
                'A' => 0b00,
                'C' => 0b01,
                'G' => 0b10,
                'T' => 0b11,
                _ => panic!("Invalid nucleotide: {}", nucleotide),
            };
            current_byte |= bits << shift; // shift the bits to the correct position

            if shift == 0 {
                data.push(current_byte);
                current_byte = 0;
                shift = 6; // reset the shift for the next byte
            } else {
                shift -= 2;
            }

            // if we are at the last nucleotide but haven't filled the current byte
            if i == dna.len() - 1 && shift != 6 {
                data.push(current_byte); // push the last byte even if partially filled
            }
        }

        DnaSequence {
            data,
            original_length: dna.len(),
        }
    }

    // decode the bit-packed DnaSequence back into a DNA string
    pub fn decode(&self) -> String {
        let mut dna = String::with_capacity(self.original_length);
        let mut shift = 6; // Start with the highest 2 bits
        let mut count = 0;

        for &byte in &self.data {
            while shift >= 0 && count < self.original_length {
                let bits = (byte >> shift) & 0b11; // Extract 2 bits
                let nucleotide = match bits {
                    0b00 => 'A',
                    0b01 => 'C',
                    0b10 => 'G',
                    0b11 => 'T',
                    _ => unreachable!(),
                };
                dna.push(nucleotide);
                count += 1;
                shift -= 2; // move to the next 2 bits
            }
            shift = 6; // reset shift for the next byte
        }

        dna
    }

    // return the original number of nucleotides
    pub fn len(&self) -> usize {
        self.original_length
    }

    pub fn get_range(&self, start: usize, end: usize) -> Vec<u8> {
        let mut nucleotides = Vec::with_capacity(end - start);

        for i in start..end {
            let byte_index = i / 4;
            let shift = 6 - (i % 4) * 2;
            let bits = (self.data[byte_index] >> shift) & 0b11;

            let nucleotide = match bits {
                0b00 => b'A',
                0b01 => b'C',
                0b10 => b'G',
                0b11 => b'T',
                _ => unreachable!(),
            };
            nucleotides.push(nucleotide);
        }
        nucleotides
    }
}
