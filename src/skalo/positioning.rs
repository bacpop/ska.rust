use hashbrown::{HashMap, HashSet};
use seq_io::fasta::{Reader, Record};
use std::cmp::Ordering;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use flate2::read::MultiGzDecoder;

use crate::skalo::utils::{encode_kmer, rev_compl, VariantInfo};

// extract genomic k-mers with up to 3 distinct positions
pub fn extract_genomic_kmers(
    file_path: PathBuf,
    k: usize,
) -> (HashMap<u128, Vec<u32>>, Vec<u8>, String) {
    // Initialize HashMap to store k-mers and their positions
    let mut kmer_map: HashMap<u128, Vec<u32>> = HashMap::new();

    // initialize HashSet to track k-mers that have more than 5 positions
    let mut overflow_kmers: HashSet<u128> = HashSet::new();

    // initialize variables to store the genome
    let mut genome_seq: Vec<u8> = Vec::new();
    let mut genome_name: String = "".to_string();

    // set the reader for compressed or uncompressed files
    let buf = get_reader(&file_path);
    let mut reader = Reader::new(buf);

    // initialize sequence counter
    let mut sequence_count = 0;

    // process records one by one
    while let Some(record) = reader.next() {
        // increment sequence counter
        sequence_count += 1;
        if sequence_count > 1 {
            panic!("\nError: more than one sequence detected in the reference genome file.\n");
        }

        // unwrap record (contains name, sequence, and quality)
        let record_ready = match record {
            Ok(record) => record,
            Err(_) => continue, // Skip on error
        };

        // process the genome sequence
        genome_seq = record_ready
            .seq()
            .iter()
            .copied()
            .filter(|&byte| !byte.is_ascii_whitespace())
            .map(|byte| byte.to_ascii_uppercase())
            .collect();

        // get the genome name
        genome_name = record_ready.id().unwrap().to_string();

        // only consider sequences long enough to have a k-mer
        if genome_seq.len() >= k {
            // extract k-mers (slices from Vec<u8>)
            for n in 0..(genome_seq.len() - k + 1) {
                // get slice of Vec<u8>
                let kmer = &genome_seq[n..n + k];

                // convert k-mer to u128
                if let Some(kmer_encoded) = encode_vecu8_u128(kmer) {
                    // Skip k-mers that have already overflowed
                    if overflow_kmers.contains(&kmer_encoded) {
                        continue;
                    }

                    // insert or update the k-mer in the HashMap
                    let positions = kmer_map.entry(kmer_encoded).or_insert_with(Vec::new);
                    if positions.len() < 3 {
                        positions.push((n + k) as u32);
                    }

                    // if positions exceed 3, remove the k-mer and add it to the overflow set
                    if positions.len() > 3 {
                        kmer_map.remove(&kmer_encoded);
                        overflow_kmers.insert(kmer_encoded);
                    }
                }
            }
        }
    }

    (kmer_map, genome_seq, genome_name)
}

// reader function
fn get_reader(path: &PathBuf) -> Box<dyn BufRead + Send> {
    let filename_str = path.to_str().unwrap();
    let file = match File::open(path) {
        Ok(file) => file,
        Err(error) => panic!("Error opening file: {:?}.", error),
    };

    if filename_str.ends_with(".gz") {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    }
}

// encode slice [u8] into a u128
fn encode_vecu8_u128(kmer: &[u8]) -> Option<u128> {
    let mut encoded = 0u128;
    for &base in kmer {
        encoded = (encoded << 2)
            | match base {
                b'A' => 0b00,
                b'C' => 0b01,
                b'G' => 0b10,
                b'T' => 0b11,
                _ => return None,
            };
    }
    Some(encoded)
}

// returns the genomic position of a bubble
pub fn scan_variants(
    vec_variants: &[VariantInfo],
    len_kmer_graph: usize,
    kmer_map: &HashMap<u128, Vec<u32>>,
) -> (bool, u32, String) {
    let mut final_position = 0;

    let mut vec_position_forward: Vec<u32> = Vec::new();
    let mut vec_position_reverse: Vec<u32> = Vec::new();

    for variant in vec_variants {
        let seq = variant.sequence.decode();
        let rc_seq = rev_compl(&seq);

        // collect positions for the forward sequence
        for pos in 0..=seq.len() - len_kmer_graph {
            let encoded_kmer = encode_kmer(&seq[pos..pos + len_kmer_graph]);
            if let Some(vec_pos) = kmer_map.get(&encoded_kmer) {
                for position in vec_pos {
                    vec_position_forward.push(position - pos as u32);
                }
            }
        }

        // collect positions for the reverse-complement sequence
        for pos in 0..=rc_seq.len() - len_kmer_graph {
            let encoded_kmer = encode_kmer(&rc_seq[pos..pos + len_kmer_graph]);
            if let Some(vec_pos) = kmer_map.get(&encoded_kmer) {
                for position in vec_pos {
                    vec_position_reverse.push(position - pos as u32);
                }
            }
        }
    }

    // determine the most frequent position and its count for forward and reverse-complement
    let forward = if !vec_position_forward.is_empty() {
        let result = most_frequent_position(&vec_position_forward);
        if result.1 == 0 {
            None
        } else {
            Some(result)
        }
    } else {
        None
    };

    let reverse = if !vec_position_reverse.is_empty() {
        let result = most_frequent_position(&vec_position_reverse);
        if result.1 == 0 {
            None
        } else {
            Some(result)
        }
    } else {
        None
    };

    // determine the final position and orientation
    let (positioned, final_orientation) = match (forward, reverse) {
        (Some((pos_forward, count_forward)), Some((pos_reverse, count_reverse))) => {
            match count_forward.cmp(&count_reverse) {
                std::cmp::Ordering::Equal => (false, "none".to_string()),
                std::cmp::Ordering::Greater => {
                    final_position = pos_forward;
                    (true, "for".to_string())
                }
                std::cmp::Ordering::Less => {
                    final_position = pos_reverse;
                    (true, "rc".to_string())
                }
            }
        }
        (Some((pos_forward, _)), None) => {
            final_position = pos_forward;
            (true, "for".to_string())
        }
        (None, Some((pos_reverse, _))) => {
            final_position = pos_reverse;
            (true, "rc".to_string())
        }
        (None, None) => (false, "none".to_string()),
    };

    (positioned, final_position, final_orientation)
}

// returns the most frequent position and its count, or (0, 0) if there's a tie
fn most_frequent_position(numbers: &[u32]) -> (u32, usize) {
    let counts = numbers.iter().fold(HashMap::new(), |mut counts, &num| {
        *counts.entry(num).or_insert(0) += 1;
        counts
    });

    let mut max_element = None;
    let mut max_count = 0;
    let mut tie = false;

    for (&num, &count) in &counts {
        match count.cmp(&max_count) {
            Ordering::Greater => {
                max_element = Some(num);
                max_count = count;
                tie = false;
            }
            Ordering::Equal => {
                tie = true;
            }
            Ordering::Less => {}
        }
    }

    if tie {
        return (0, 0); // return (0, 0) in case of a tie
    }

    // if no tie, return the most frequent position and its count
    if let Some(position) = max_element {
        if max_count < 10 {
            return (0, 0); // return (0, 0) in case of low positioning
        } else {
            return (position, max_count);
        }
    }

    (0, 0) // default return if input is empty
}
