use hashbrown::{HashMap, HashSet};
use std::fs::File;
use std::io::Write;

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

pub fn identify_good_kmers(
    len_kmer: usize,
    all_kmers: &HashMap<u128, HashMap<u128, u32>>,
    index_map: &HashMap<u32, String>,
) -> (HashSet<u128>, HashSet<u128>) {
    println!(" # identify bubble extremities");

    let len_kmer_graph = len_kmer - 1;

    let mut start_kmers: HashSet<u128> = HashSet::new();
    let mut end_kmers: HashSet<u128> = HashSet::new();

    for (kmer, next_kmers_map) in all_kmers.iter() {
        if next_kmers_map.len() > 1 {
            let all_next_kmer: Vec<&u128> = next_kmers_map.keys().collect();

            'i_loop: for (i, &next_kmer) in all_next_kmer.iter().enumerate() {
                for &next_kmer2 in all_next_kmer.iter().skip(i + 1) {
                    let samples_1 = &index_map[&next_kmers_map[next_kmer]];
                    let samples_2 = &index_map[&next_kmers_map[next_kmer2]];

                    if compare_samples(samples_1, samples_2) {
                        start_kmers.insert(kmer.clone());

                        let dna = decode_kmer(kmer.clone(), len_kmer_graph);
                        let rc = rev_compl(&dna);

                        //uncomment to print network
                        //println!("{}	{}	red", &dna, &dna);

                        end_kmers.insert(encode_kmer(&rc));

                        //uncomment to print network
                        //println!("{}	{}	red", &rc, &rc);
                        break 'i_loop;
                    }
                }
            }
        }
    }

    // exit program if no extremity found (e.g. cases of weeded skf files)
    if start_kmers.is_empty() {
        eprintln!("\n      Error: there is no entry node in this graph, hence no variant.\n");
        std::process::exit(1);
    }

    println!("     . {} entry nodes", start_kmers.len());
    (start_kmers, end_kmers)
}

fn compare_samples(str_1: &str, str_2: &str) -> bool {
    let set1: HashSet<&str> = str_1.split('|').collect();
    let set2: HashSet<&str> = str_2.split('|').collect();
    let diff_1 = set1.difference(&set2).count();
    let diff_2 = set2.difference(&set1).count();

    diff_1 != 0 && diff_2 != 0
}

pub fn filter_output_sequences(
    variant_groups: HashMap<String, Vec<VariantInfo>>,
    len_kmer: usize,
    all_samples: Vec<String>,
    n_homopolymer: Option<u32>,
    max_missing: f32,
    output_name: &str,
    input_name: &str,
) {
    println!(" # filter sequences");

    let len_kmer_graph = len_kmer - 1;

    // prepare output files
    let mut out_1 =
        File::create(format!("{}_seq_groups.fas", output_name)).expect("Failed to create file");
    let mut out_2 =
        File::create(format!("{}_variants.tsv", output_name)).expect("Failed to create file");

    // write header of TSV variant file
    out_2
        .write_all(format!("#skalo version {}\n", env!("CARGO_PKG_VERSION")).as_bytes())
        .expect("Failed to write to file");
    out_2
        .write_all(
            format!(
                "#parameters: m={}{}\n",
                max_missing,
                match n_homopolymer {
                    Some(n) => format!(" n={}", n),
                    None => String::new(),
                }
            )
            .as_bytes(),
        )
        .expect("Failed to write to file");
    out_2
        .write_all(format!("#input file: {} (k={})\n", input_name, len_kmer).as_bytes())
        .expect("Failed to write to file");

    let index_name: Vec<String> = all_samples
        .iter()
        .enumerate()
        .map(|(index, value)| format!("{}:{}", index, value))
        .collect();
    let formatted_string = index_name.join(", ");
    out_2
        .write_all(format!("#samples: {}\n", formatted_string).as_bytes())
        .expect("Failed to write to file");
    out_2.write_all(b"#pos_ali\tnb_states\ttype\tunclear_insert\tfirst_kmer\tvariants\tlast_kmer\tratio_missing\tsamples\n")
        .expect("Failed to write to file");

    // Prepare variables for binary alignment
    let mut binary_seq: HashMap<String, Vec<String>> = HashMap::new();
    let mut d_samples: HashMap<String, String> = HashMap::new();
    for (i, sample) in all_samples.iter().enumerate() {
        binary_seq.insert(sample.clone(), Vec::new());
        d_samples.insert(i.to_string(), sample.clone()); // map sample ID to sample full name
    }

    let mut nb_missing = 0;
    let mut nb_homopolymer = 0;

    // variant groups 1 by 1
    let mut position = 0;
    for vec_variants in variant_groups.values() {
        if vec_variants[0].is_snp {

            // do something with SNPs
        } else {
            // check the ratio of missing samples
            let ratio_missing = calculate_ratio_missing(all_samples.len(), vec_variants.clone());

            if ratio_missing > max_missing {
                nb_missing += 1;
            } else {
                // CHECK IF UNCLEAR INSERT
                // extract 1st kmer using 1st sequence
                let first_kmer = vec_variants[0].sequence[..len_kmer_graph].to_string();

                // collect 'middle-bases' and last k-mer (forward and reverse-complement)
                let (l_middle_bases, last_kmer) =
                    collect_middle_bases(vec_variants.clone(), len_kmer_graph, false);
                let (rc_l_middle_bases, rc_last_kmer) =
                    collect_middle_bases(vec_variants.clone(), len_kmer_graph, true);

                // test if multiple positions exist for the last k-mer
                let multiple_positions =
                    test_multiple_positions(l_middle_bases.clone(), last_kmer.clone());
                let multiple_positions_rc =
                    test_multiple_positions(rc_l_middle_bases, rc_last_kmer);

                let mut unclear_insert = "no";
                if multiple_positions || multiple_positions_rc {
                    unclear_insert = "yes";
                }

                // CHECK IK HOMOPOLYMER
                // check if indel in homopolymer
                let mut is_homopolymer = false;
                if let Some(mut n_max) = n_homopolymer {
                    // adjust n_homopolymer if higher than k-mer length to avoid the program crashing
                    if n_max > len_kmer_graph.try_into().unwrap() {
                        n_max = len_kmer_graph as u32;
                        println!("     . n_homopolymer was reduced to fit the k-mer length");
                    }
                    // test presence homopolymer
                    is_homopolymer =
                        check_homopolymer(n_max, first_kmer.clone(), l_middle_bases.clone());
                }

                if is_homopolymer {
                    nb_homopolymer += 1;
                } else {
                    // Output seq_groups
                    for variant in vec_variants.iter() {
                        let str_samples = variant
                            .maj_samples
                            .iter()
                            .copied()
                            .collect::<Vec<_>>()
                            .join("|");
                        out_1
                            .write_all(
                                format!(">{}_{}\n{}\n", position, str_samples, variant.sequence)
                                    .as_bytes(),
                            )
                            .expect("Failed to write to file");
                    }

                    // update binary alignment ('-' = missing data; '?' = both states (= 'unknown'))
                    // only consider the 2 most frequent variants in cases of 3+
                    let mut sample_done: HashSet<String> = HashSet::new();
                    let mut state = 0;
                    for variant in vec_variants.iter().take(2) {
                        // collect sample names for this sequence and update its vector in binary_seq
                        for sample_id in &variant.maj_samples {
                            let full_name = d_samples.get(sample_id.to_owned()).unwrap();
                            sample_done.insert(full_name.clone());
                            let sample_vec = binary_seq.get_mut(full_name).unwrap();

                            if sample_vec.len() <= position {
                                sample_vec.resize(position + 1, state.to_string());
                            } else {
                                sample_vec[position] = "?".to_string();
                            }
                        }
                        // update number of character states
                        state += 1;
                    }

                    // update binary alignment with missing samples
                    for sample in &all_samples {
                        if !sample_done.contains(sample) {
                            binary_seq
                                .entry(sample.to_string())
                                .and_modify(|vec| vec.push("-".to_string()));
                        }
                    }

                    // get type of variant
                    let mut type_variant = "complex";
                    for xx in &l_middle_bases {
                        if xx == "." {
                            type_variant = "_indel_";
                        }
                    }

                    // get lists of samples
                    let mut l_samples = Vec::new();
                    for variant in vec_variants.iter() {
                        // sort list of samples by equivalent integers
                        let mut vect_idx = variant.maj_samples.iter().copied().collect::<Vec<_>>();
                        vect_idx.sort_by(|a, b| {
                            let a_int = a.parse::<i32>().unwrap_or(i32::MAX);
                            let b_int = b.parse::<i32>().unwrap_or(i32::MAX);
                            a_int.cmp(&b_int)
                        });
                        // save it
                        l_samples.push(vect_idx.join(","));
                    }

                    // round ratio_missing to the 2nd decimal
                    let rounded_missing = (ratio_missing * 100.0).round() / 100.0;

                    // save variants information
                    out_2
                        .write_all(
                            format!(
                                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                                position,
                                l_middle_bases.len(),
                                type_variant,
                                unclear_insert,
                                first_kmer,
                                l_middle_bases.join(" / "),
                                last_kmer,
                                rounded_missing,
                                l_samples.join(" / ")
                            )
                            .as_bytes(),
                        )
                        .expect("Failed to write to file");

                    // Update position
                    position += 1;
                }
            }
        }
    }

    // Close output files
    out_1.flush().expect("Failed to flush file");
    out_2.flush().expect("Failed to flush file");

    // Output alignment
    let mut out_3 =
        File::create(format!("{}_binary_ali.fas", output_name)).expect("Failed to create file");
    for (sample, l_seq) in &binary_seq {
        out_3
            .write_all(format!(">{}\n{}\n", sample, l_seq.join("")).as_bytes())
            .expect("Failed to write to file");
    }
    out_3.flush().expect("Failed to flush file");

    // final message
    println!("     . {} removed because of missing data", nb_missing);
    if n_homopolymer.is_some() {
        println!("     . {} removed because in homopolymers", nb_homopolymer);
    }
    println!("     . {} FINAL variant groups", position);
    println!("done.");
}

fn calculate_ratio_missing(nb_total: usize, vec_variants: Vec<VariantInfo>) -> f32 {
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

fn collect_middle_bases(
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

fn test_multiple_positions(l_middles: Vec<String>, last_kmer: String) -> bool {
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

fn check_homopolymer(limit_n: u32, first_kmer: String, l_middle_bases: Vec<String>) -> bool {
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
