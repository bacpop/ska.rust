use bit_set::BitSet;
use hashbrown::{HashMap, HashSet};

use crate::ska_dict::bit_encoding::UInt;
use crate::skalo::output::create_fasta_and_vcf;
use crate::skalo::positioning::{extract_genomic_kmers, scan_variants};
use crate::skalo::utils::{Config, DataInfo, VariantInfo};

type VariantGroups<IntT> = HashMap<(IntT, IntT), Vec<VariantInfo>>;

pub fn analyse_variant_groups<IntT: for<'a> UInt<'a>>(
    mut variant_groups: VariantGroups<IntT>,
    indel_groups: VariantGroups<IntT>,
    kmer_2_samples: HashMap<IntT, BitSet>,
    config: &Config,
    data_info: &DataInfo,
) {
    // check if the optional reference genome file argument is provided -> extract kmers
    let (do_postioning, kmer_map, genome_name, genome_seq) =
        if let Some(path) = &config.reference_genome {
            log::info!("Reading reference genome");
            let (extracted_kmer_map, seq, name) =
                extract_genomic_kmers(path.clone(), data_info.k_graph);
            (true, extracted_kmer_map, name, seq)
        } else {
            (
                false,
                HashMap::<u128, Vec<u32>>::new(),
                "".to_string(),
                Vec::<u8>::new(),
            )
        };

    log::info!("Processing indels");

    // collect entry kmers of indels
    let (final_indels, entries_indels) = process_indels(indel_groups, data_info.k_graph);
    log::info!("Found {} indels", final_indels.len());

    // remove variants having  internal indels from each variant group
    for (_, vec_variant) in variant_groups.iter_mut() {
        let mut i = 0;
        while i < vec_variant.len() {
            let nb_indel_kmers = find_internal_indels(&vec_variant[i], &entries_indels, data_info);
            // there has to be 4 ends for 2 indels, but reducing the threshold to 3 half the numbers of FPs
            if nb_indel_kmers > config.max_indel_kmers {
                vec_variant.remove(i);
            } else {
                i += 1;
            }
        }
    }

    log::info!("Processing SNPs");

    // create a vector of keys sorted by the ratio of size of Vec<VariantInfo> to the length of the first sequence
    // and sort the keys by decreasing order -> we consider first for snp calling variant group with lot of variants
    let mut sorted_keys: Vec<_> = variant_groups
        .iter()
        .filter_map(|(key, value)| {
            // Access the first sequence length safely
            value.first().map(|variant_info| {
                let sequence_length = variant_info.sequence.len();
                let ratio = value.len() as f64 / sequence_length as f64;
                (key, ratio)
            })
        })
        .collect();
    sorted_keys.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap()); // Sort by ratio, descending

    // start processing SNPs
    let mut entries_done: HashSet<IntT> = HashSet::new();

    // to store SNPs, with genomic position as key (or counter if no positioning)
    let mut final_snps: HashMap<u32, Vec<char>> = HashMap::new();
    let mut not_postioned = 0;
    let mut counter = 0;

    for (key, _) in sorted_keys {
        if !entries_indels.contains(&key.0)
            && !entries_indels.contains(&IntT::rev_comp(key.1, data_info.k_graph))
        {
            let vec_variants = variant_groups.get(key).unwrap();

            // case with 0 or 1 seq left (ie, not enough to be a variant group)
            if vec_variants.len() < 2 {
                continue;
            }

            // check potential SNPs from VariantInfo
            let real_snp_pos = get_potential_snp(vec_variants);

            // get SNP column and kmers
            let mut kmers_to_save: HashSet<IntT> = HashSet::new();
            let mut found_snp_pos: HashMap<usize, Vec<char>> =
                HashMap::with_capacity(real_snp_pos.len());

            for &pos in &real_snp_pos {
                let mut snp_column = vec!['-'; data_info.sample_names.len()];
                let mut tmp_kmers: HashSet<IntT> = HashSet::new();

                let mut new_snp = true;

                for variant in vec_variants {
                    let seq = &variant.sequence;

                    // Extract k-mers directly from the packed DNA sequence
                    let full_before =
                        IntT::encode_kmer(&seq.get_range(pos - data_info.k_graph, pos + 1));
                    let full_after =
                        IntT::encode_kmer(&seq.get_range(pos, pos + data_info.k_graph + 1));
                    let rc_after = IntT::rev_comp(full_after, data_info.k_graph + 1);

                    // this is the critical part: we have to avoid SNPs already identified
                    if !entries_done.contains(&full_before) && !entries_done.contains(&rc_after) {
                        let last_nucl = IntT::get_last_nucl(full_before);
                        let samples = kmer_2_samples.get(&full_before).unwrap();

                        for sample_index in samples {
                            if snp_column[sample_index] == '-'
                                || snp_column[sample_index] == last_nucl
                            {
                                snp_column[sample_index] = last_nucl;
                            } else {
                                snp_column[sample_index] = 'N';
                            }
                        }

                        // Save k-mers to avoid
                        tmp_kmers.insert(full_before);
                        tmp_kmers.insert(IntT::rev_comp(full_before, data_info.k_graph + 1));
                        tmp_kmers.insert(full_after);
                        tmp_kmers.insert(rc_after);
                    } else {
                        new_snp = false;
                    }
                }
                // check level of missing data if new SNP
                if new_snp {
                    let (true_variant, ratio_missing) =
                        check_missing_data(data_info.sample_names.len(), &snp_column);
                    if true_variant && ratio_missing <= config.max_missing {
                        // save surrounding k-mers
                        kmers_to_save.extend(tmp_kmers);
                        // save SNP
                        found_snp_pos.insert(pos, snp_column);
                    }
                }
            }
            entries_done.extend(kmers_to_save.iter());

            // variant positioning if reference genome and if a SNP has been found
            if !found_snp_pos.is_empty() {
                if do_postioning {
                    let (position_found, position, orientation) =
                        scan_variants(vec_variants, data_info.k_graph, &kmer_map);

                    if position_found {
                        let seq_length = vec_variants[0].sequence.len();
                        let is_forward = orientation == "for";

                        // adjust position with SNP pos in variant group and orientation
                        for (pos, column) in found_snp_pos {
                            let final_position = if is_forward {
                                position + (pos - data_info.k_graph) as u32
                            } else {
                                position + (seq_length - pos - data_info.k_graph - 1) as u32
                            };

                            let final_column = if is_forward {
                                column
                            } else {
                                complement_snp(&column)
                            };

                            // save it if position not already taken
                            if final_snps.contains_key(&final_position) {
                                not_postioned += 1;
                            } else {
                                final_snps.insert(final_position, final_column);
                            }
                        }
                    } else {
                        not_postioned += found_snp_pos.len();
                    }
                } else {
                    for (_, column) in found_snp_pos {
                        counter += 1;
                        // save it
                        final_snps.insert(counter, column);
                    }
                }
            }
        }
    }

    if do_postioning {
        log::info!(
            "     -> {} SNPs (+ {} w/o position)",
            final_snps.len(),
            not_postioned
        );
    } else {
        log::info!("     -> {} SNPs", final_snps.len());
    }

    log::info!("Writing output");
    // write output
    create_fasta_and_vcf(
        genome_name,
        genome_seq,
        data_info.sample_names.clone(),
        final_snps,
        config,
    );
}

fn find_internal_indels<IntT: for<'a> UInt<'a>>(
    variant: &VariantInfo,
    entries_indels: &HashSet<IntT>,
    data_info: &DataInfo,
) -> usize {
    let mut nb = 0;
    let sequence = &variant.sequence;
    let k_graph = data_info.k_graph;

    // precompute the initial k-mer
    let mut kmer = IntT::encode_kmer(&sequence.get_range(0, k_graph));
    // Mask for retaining k-mer length
    let mask = IntT::skalo_mask(k_graph);

    // sliding window for k-mer computation
    for i in k_graph..sequence.len() {
        // compute the next k-mer using a rolling hash
        let byte_index = i / 4;
        let shift = 6 - (i % 4) * 2;
        let next_nucleotide = IntT::from_encoded_base((sequence.data[byte_index] >> shift) & 0b11);

        kmer = ((kmer << 2) | next_nucleotide) & mask; // update k-mer with new nucl

        if entries_indels.contains(&kmer) {
            nb += 1;
        }
    }
    nb
}

fn process_indels<IntT: for<'a> UInt<'a>>(
    indel_groups: VariantGroups<IntT>,
    k_graph: usize,
) -> (VariantGroups<IntT>, HashSet<IntT>) {
    let mut entries_indels: HashSet<IntT> = HashSet::new();
    let mut final_indels: VariantGroups<IntT> = HashMap::new();

    for (combined_ext, vec_variant) in &indel_groups {
        if !entries_indels.contains(&combined_ext.0) {
            // test if rev-compl exists
            let rc_1 = IntT::rev_comp(combined_ext.0, k_graph);
            let rc_2 = IntT::rev_comp(combined_ext.1, k_graph);
            let rc_combined = (rc_2, rc_1);

            if indel_groups.contains_key(&rc_combined) {
                // compare the size of variant groups and select the shortest one
                // this is equivalent to indel realigning in read-alignment (useful in repeats)
                let vec_variant2 = indel_groups.get(&rc_combined).unwrap();
                let sum1: usize = vec_variant
                    .iter()
                    .map(|variant| variant.sequence.len())
                    .sum();
                let sum2: usize = vec_variant2
                    .iter()
                    .map(|variant| variant.sequence.len())
                    .sum();
                if sum1 <= sum2 {
                    final_indels.insert(*combined_ext, vec_variant.clone());
                } else {
                    final_indels.insert(rc_combined, vec_variant2.clone());
                }
            } else {
                final_indels.insert(*combined_ext, vec_variant.clone());
            }
            entries_indels.insert(combined_ext.0);
            entries_indels.insert(rc_1);
            entries_indels.insert(combined_ext.1);
            entries_indels.insert(rc_2);
        }
    }
    (final_indels, entries_indels)
}

fn get_potential_snp(vec_variant: &Vec<VariantInfo>) -> HashSet<usize> {
    let mut snps_set = HashSet::new();
    // collect all SNPs into the HashSet
    for variant in vec_variant {
        snps_set.extend(&variant.vec_snps);
    }

    let mut actual_snps_set = HashSet::new();

    // check which positions in snps_set are actual SNPs
    for &pos in &snps_set {
        let mut nucleotide_presence = [false; 4]; // A, C, G, T

        for variant in vec_variant {
            let seq = &variant.sequence;
            if pos < seq.len() {
                let nucl = seq.get_range(pos, pos + 1)[0];
                match nucl {
                    b'A' => nucleotide_presence[0] = true,
                    b'C' => nucleotide_presence[1] = true,
                    b'G' => nucleotide_presence[2] = true,
                    b'T' => nucleotide_presence[3] = true,
                    _ => {}
                }
            }
        }

        // count the number of distinct nucleotides
        let distinct_count = nucleotide_presence.iter().filter(|&&x| x).count();
        if distinct_count > 1 {
            actual_snps_set.insert(pos);
        }
    }
    actual_snps_set
}

pub fn check_missing_data(nb_total: usize, snp_column: &[char]) -> (bool, f32) {
    // count occurrences of valid SNPs (A, T, G, C) and calculate missing data
    let mut nucleotide_counts = [false; 4];
    let mut missing_samples = 0;

    for &snp in snp_column {
        match snp {
            'A' => nucleotide_counts[0] = true,
            'T' => nucleotide_counts[1] = true,
            'G' => nucleotide_counts[2] = true,
            'C' => nucleotide_counts[3] = true,
            _ => missing_samples += 1,
        }
    }

    // calculate ratio of missing samples
    let ratio_missing = missing_samples as f32 / nb_total as f32;

    // check for at least two nucleotides present among A, T, G, C
    let valid_nucleotide_count = nucleotide_counts.iter().filter(|&&present| present).count();

    (valid_nucleotide_count >= 2, ratio_missing)
}

fn complement_snp(dna: &[char]) -> Vec<char> {
    dna.iter()
        .map(|&nucleotide| match nucleotide {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            '-' => '-',
            'N' => 'N',
            _ => panic!("Invalid nucleotide: {}", nucleotide),
        })
        .collect()
}
