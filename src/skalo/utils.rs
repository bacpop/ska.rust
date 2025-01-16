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
    pub	max_indel_kmers: usize,
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




/*
pub fn resolve_upac(nucl1: &str, nucl2: &str) -> String {
    match (nucl1, nucl2) {
        // Basic nucleotide conflicts
        ("A", "G") | ("G", "A") => "R".to_string(),
        ("C", "T") | ("T", "C") => "Y".to_string(),
        ("G", "C") | ("C", "G") => "S".to_string(),
        ("A", "T") | ("T", "A") => "W".to_string(),
        ("G", "T") | ("T", "G") => "K".to_string(),
        ("A", "C") | ("C", "A") => "M".to_string(),

        // Handling conflicts with existing UPAC codes
        ("A", "R") | ("R", "A") | ("G", "R") | ("R", "G") => "R".to_string(),
        ("C", "Y") | ("Y", "C") | ("T", "Y") | ("Y", "T") => "Y".to_string(),
        ("G", "S") | ("S", "G") | ("C", "S") | ("S", "C") => "S".to_string(),
        ("A", "W") | ("W", "A") | ("T", "W") | ("W", "T") => "W".to_string(),
        ("G", "K") | ("K", "G") | ("T", "K") | ("K", "T") => "K".to_string(),
        ("A", "M") | ("M", "A") | ("C", "M") | ("M", "C") => "M".to_string(),

        // Mixed UPAC code resolution with specific nucleotides
        ("A", "Y") | ("Y", "A") => "V".to_string(),
        ("G", "Y") | ("Y", "G") => "B".to_string(),
        ("T", "R") | ("R", "T") => "D".to_string(),
        ("C", "R") | ("R", "C") => "V".to_string(),

        // Extended UPAC combinations for three nucleotides
        ("A", "V") | ("V", "A") | ("C", "V") | ("V", "C") | ("G", "V") | ("V", "G") => "V".to_string(),
        ("A", "H") | ("H", "A") | ("C", "H") | ("H", "C") | ("T", "H") | ("H", "T") => "H".to_string(),
        ("A", "D") | ("D", "A") | ("G", "D") | ("D", "G") | ("T", "D") | ("D", "T") => "D".to_string(),
        ("C", "B") | ("B", "C") | ("G", "B") | ("B", "G") | ("T", "B") | ("B", "T") => "B".to_string(),

        // Default to N for any unresolved conflicts
        _ => "N".to_string(),
    }
}
*/


pub fn rev_compl(seq: &str) -> String {
    let reverse_code: HashMap<char, char> = [
        ('A', 'T'),
        ('C', 'G'),
        ('G', 'C'),
        ('T', 'A'),
        ('N', 'N'),
    ]
    .iter()
    .cloned()
    .collect();

    let reversed_seq: String = seq.chars().rev().collect();
    let complemented_seq: String = reversed_seq.chars().map(|n| reverse_code[&n]).collect();

    complemented_seq
}


pub fn rev_compl_u128(kmer: u128, k: usize) -> u128 {
    // mask for the last 2 bits (representing one nucleotide)
    let mask = 0b11u128;

    // initialize reverse complement result
    let mut rc_u128 = 0u128;

    // for each nucleotide in the k-mer
    for i in 0..k {
        // get the last 2 bits (current nucleotide)
        let nucleotide = (kmer >> (2 * i)) & mask;

        // complement the nucleotide:
        let complement = match nucleotide {
            0b00 => 0b11, // A -> T
            0b01 => 0b10, // C -> G
            0b10 => 0b01, // G -> C
            0b11 => 0b00, // T -> A
            _ => unreachable!(),
        };

        // shift the complemented nucleotide to its reverse position
        rc_u128 |= complement << (2 * (k - i - 1));
    }
    rc_u128
}


pub fn encode_kmer(kmer: &str) -> u128 {
    let nucleotide_to_bits: [u8; 4] = [
        0b00, // A
        0b01, // C
        0b10, // G
        0b11, // T
    ];

    let mut result: u128 = 0;

    for nucleotide in kmer.chars() {
        let index = match nucleotide {
            'A' => 0,
            'C' => 1,
            'G' => 2,
            'T' => 3,
            _ => panic!("Invalid nucleotide"),
        };
        result = (result << 2) | (nucleotide_to_bits[index] as u128);
    }
    result
}


pub fn encode_u8_kmer(dna: &[u8]) -> u128 {
    let mut encoded: u128 = 0;
    for &nucleotide in dna {
        encoded <<= 2; // Shift left by 2 bits
        encoded |= match nucleotide {
            b'A' => 0b00,
            b'C' => 0b01,
            b'G' => 0b10,
            b'T' => 0b11,
            _ => panic!("Invalid nucleotide: {}", nucleotide as char),
        };
    }
    encoded
}

pub fn decode_kmer(encoded: u128, k: usize) -> String {
    let bits_to_nucleotide: [char; 4] = ['A', 'C', 'G', 'T'];
    let mut kmer = String::with_capacity(k);

    let mask = (1u128 << (2 * k)) - 1;
    let mut value = encoded & mask;

    for _ in 0..k {
        let index = (value & 0b11) as usize;
        let nucleotide = bits_to_nucleotide[index];
        kmer.insert(0, nucleotide);
        value >>= 2;
    }
    kmer
}


// extract last nucleotide from an encoded k-mer
pub fn get_last_nucl(encoded_kmer: u128) -> char {
    // mask the last 2 bits to get the encoded nucleotide
    let last_bits = (encoded_kmer & 0b11) as u8;
    // decode the nucleotide based on the 2-bit pattern
    match last_bits {
        0b00 => 'A',
        0b01 => 'C',
        0b10 => 'G',
        0b11 => 'T',
        _ => unreachable!(),
    }
}


#[derive(Clone)]
pub struct VariantInfo {
    pub sequence: DnaSequence,
    pub vec_snps: Vec<usize>,
}

impl VariantInfo {
    pub fn new(
        sequence: DnaSequence,
        vec_snps: Vec<usize>,
    ) -> Self {
        VariantInfo {
            sequence,
            vec_snps,
        }
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
        let mut shift = 6;  // start with the highest 2 bits

        for (i, nucleotide) in dna.chars().enumerate() {
            let bits = match nucleotide {
                'A' => 0b00,
                'C' => 0b01,
                'G' => 0b10,
                'T' => 0b11,
                _ => panic!("Invalid nucleotide: {}", nucleotide),
            };
            current_byte |= bits << shift;  // shift the bits to the correct position

            if shift == 0 {
                data.push(current_byte);
                current_byte = 0;
                shift = 6;  // reset the shift for the next byte
            } else {
                shift -= 2;
            }

            // if we are at the last nucleotide but haven't filled the current byte
            if i == dna.len() - 1 && shift != 6 {
                data.push(current_byte);  // push the last byte even if partially filled
            }
        }

        DnaSequence { data, original_length: dna.len() }
    }

    // decode the bit-packed DnaSequence back into a DNA string
    pub fn decode(&self) -> String {
        let mut dna = String::with_capacity(self.original_length);
        let mut shift = 6;  // Start with the highest 2 bits
        let mut count = 0;

        for &byte in &self.data {
            while shift >= 0 && count < self.original_length {
                let bits = (byte >> shift) & 0b11;  // Extract 2 bits
                let nucleotide = match bits {
                    0b00 => 'A',
                    0b01 => 'C',
                    0b10 => 'G',
                    0b11 => 'T',
                    _ => unreachable!(),
                };
                dna.push(nucleotide);
                count += 1;
                shift -= 2;  // move to the next 2 bits
            }
            shift = 6;  // reset shift for the next byte
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

