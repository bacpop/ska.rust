use hashbrown::{HashMap, HashSet};
use std::str::FromStr;

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

pub fn decode_kmer_skalo(encoded: u128, k: usize) -> String {
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

pub fn get_last_nucleotide(encoded_kmer: u128) -> char {
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

/// structure to save variant information

#[derive(Clone)]
pub struct VariantInfo<'a> {
    //pub entry_kmer: u128,
    //pub exit_kmer: u128,
    pub sequence: String,
    pub is_snp: bool,
    pub visited_nodes: HashSet<u128>,
    pub count_samples: HashMap<&'a str, i32>,
    pub maj_samples: HashSet<&'a str>,
}

impl<'a> VariantInfo<'a> {
    pub fn new(
        //entry_kmer: u128,
        //exit_kmer: u128,
        sequence: String,
        is_snp: bool,
        visited_nodes: HashSet<u128>,
        count_samples: HashMap<&'a str, i32>,
        maj_samples: HashSet<&'a str>,
    ) -> Self {
        VariantInfo {
            //entry_kmer,
            //exit_kmer,
            sequence,
            is_snp,
            visited_nodes,
            count_samples,
            maj_samples,
        }
    }
}

pub fn compare_samples(str_1: &str, str_2: &str) -> bool {
    let set1: HashSet<&str> = str_1.split('|').collect();
    let set2: HashSet<&str> = str_2.split('|').collect();
    let diff_1 = set1.difference(&set2).count();
    let diff_2 = set2.difference(&set1).count();

    diff_1 != 0 && diff_2 != 0
}

pub fn calculate_ratio_missing(nb_total: usize, vec_variants: Vec<VariantInfo>) -> f32 {
    // initialise vector of 0
    let mut ref_samples: Vec<u32> = vec![0; nb_total];

    // iterate through vec_variants to count present samples for all variants
    for variant in vec_variants {
        for sample in &variant.maj_samples {
            let idx = sample.parse::<usize>().unwrap();
            ref_samples[idx] += 1;
        }
    }

    // count 'missing' samples (i.e., count != 1)
    let missing_samples = ref_samples.iter().filter(|&&count| count != 1).count() as f32;

    // calculate ratio of missing samples
    let ratio_missing = missing_samples / nb_total as f32;
    ratio_missing
}

pub fn collect_middle_bases(
    vec_variants: Vec<VariantInfo>,
    len_k_graph: usize,
    use_rev_complement: bool,
) -> (Vec<String>, String) {
    // collect all sequences without the first kmer
    //let reduced_seq: Vec<&str> = vec_variants.iter().map(|variant| &variant.sequence[len_k_graph..]).collect();
    let reduced_seq: Vec<String> = vec_variants
        .iter()
        .map(|variant| {
            let sequence = if use_rev_complement {
                rev_compl(&variant.sequence)
            } else {
                variant.sequence.clone()
            };
            sequence[len_k_graph..].to_string()
        })
        .collect();

    // get start position of last k-mer (i.e. find last position for which sequences differ (from the end))
    let mut identical = true;
    let mut n_nucl = 0;

    while identical {
        n_nucl += 1;
        let mut all_ends: HashSet<String> = HashSet::new();
        // extract last n nucleotide from each seq
        for seq in &reduced_seq {
            if n_nucl > seq.len() {
                identical = false;
            } else {
                let last_n_chars: Vec<String> = seq
                    .chars()
                    .rev()
                    .take(n_nucl)
                    .map(|c| c.to_string())
                    .collect();
                let concatenated_last_chars: String = last_n_chars.into_iter().rev().collect();
                all_ends.insert(concatenated_last_chars.clone());
            }
        }
        if all_ends.len() > 1 {
            identical = false;
        }
    }
    n_nucl -= 1;

    // extract last kmer (using first sequence)
    let pos_end = reduced_seq[0].len() - n_nucl;
    let mut last_kmer = reduced_seq[0][pos_end..].to_string();

    // the length of last k-mer might be in some very rare cases longer than expected (only observed in variants with lot of missing samples) -> truncate it
    if last_kmer.len() > len_k_graph {
        last_kmer = last_kmer[..len_k_graph].to_string();
    }

    // extract 'middle-bases' (remove last kmer from reduced sequences -> only middle base left)
    let mut l_middles: Vec<String> = Vec::new();
    for seq in &reduced_seq {
        let pos2_end = seq.len() - n_nucl;
        let mut middle_bases = &seq[..pos2_end];
        if middle_bases == "" {
            middle_bases = ".";
        }
        l_middles.push(middle_bases.to_string());
    }

    (l_middles, last_kmer)
}

pub fn test_multiple_positions(l_middles: Vec<String>, last_kmer: String) -> bool {
    for middle in l_middles {
        let full_string = middle.to_string() + &last_kmer;

        let mut occurrences_count = 0;
        let mut start = 0;

        while let Some(index) = full_string[start..].find(&last_kmer) {
            let full_index = start + index;
            occurrences_count += 1;
            if occurrences_count > 1 {
                return true;
            }
            start = full_index + 1;
        }
    }

    false // only one occurrence found
}

pub fn check_homopolymer(limit_n: u32, first_kmer: String, l_middle_bases: Vec<String>) -> bool {
    let mut is_homopolymer = false;

    // extract last limit_n nucleotides of first k-mer
    let last_n_chars = &first_kmer[first_kmer.len() - limit_n as usize..];

    // check if these last limit_n nucleotides are the same (ie, homopolymer)
    let char_set: HashSet<_> = last_n_chars.chars().collect();
    if char_set.len() == 1 {
        // check if the insert correspond to the homopolymer repeat
        'loop_middle: for middle in l_middle_bases {
            let middle_char_set: HashSet<_> = middle.chars().collect();
            if middle_char_set == char_set {
                is_homopolymer = true;
                break 'loop_middle;
            }
        }
    }
    is_homopolymer
}

pub fn jaccard_similarity(set1: &HashSet<&str>, set2: &HashSet<&str>) -> f32 {
    // returns the Jaccard similarity value between 2 sets
    let intersection_size = set1.intersection(set2).count() as f32;
    let union_size = (set1.len() as f32 + set2.len() as f32 - intersection_size) as f32;
    intersection_size / union_size
}

pub fn count_occurrences(vec: &Vec<u32>) -> HashMap<u32, i32> {
    let mut count_map = HashMap::new();
    for num in vec {
        *count_map.entry(*num).or_insert(0) += 1;
    }
    count_map
}

pub fn convert_combined(input: &str, k: usize) -> Result<String, String> {
    let parts: Vec<&str> = input.split('@').collect();

    // parse the u128 values from the string parts
    let u128_1 = u128::from_str(parts[0]).map_err(|e| e.to_string())?;
    let u128_2 = u128::from_str(parts[1]).map_err(|e| e.to_string())?;

    // compute the reverse complements for both u128 values using the same k
    let rev_compl_1 = rev_compl_u128(u128_1, k);
    let rev_compl_2 = rev_compl_u128(u128_2, k);

    Ok(format!("{}@{}", rev_compl_2, rev_compl_1))
}
