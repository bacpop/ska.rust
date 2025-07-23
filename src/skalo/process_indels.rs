//! indel processing
use bit_set::BitSet;
use hashbrown::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Write};

use crate::ska_dict::bit_encoding::UInt;
use crate::skalo::utils::{Config, DataInfo, VariantInfo};

type VariantGroups<IntT> = HashMap<(IntT, IntT), Vec<VariantInfo>>;

/// This function processes indels. These are dereplicated and inserts are extracted.
/// As for SNPs, indels are filtered to only retain true variants and based on the
/// proportion of missing samples.
pub fn process_indels<IntT: for<'a> UInt<'a>>(
    indel_groups: VariantGroups<IntT>,
    kmer_2_samples: &HashMap<IntT, BitSet>,
    data_info: &DataInfo,
    config: &Config,
) -> HashSet<IntT> {
    log::info!("Processing indels");

    // dereplicate indels
    let (final_indels, entries_indels) = dereplicate_indels(indel_groups, data_info.k_graph);

    // create VCF output file
    let vcf_filename = format!("{}_indels.vcf", config.output_name);
    let file = File::create(&vcf_filename).expect("Unable to create VCF file");
    let mut writer = BufWriter::new(file);

    // Wwrite VCF header
    writeln!(writer, "##fileformat=VCFv4.2").unwrap();
    writeln!(
        writer,
        "# REF corresponds to the most frequent variant among samples"
    )
    .unwrap();
    writeln!(
        writer,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}",
        data_info.sample_names.join("\t")
    )
    .unwrap();

    let mut nb_indels = 0;

    // consider indels 1 by one
    for vec_variants in final_indels.values() {
        // get taxonomic sampling for each variant
        let bitset_vec: Vec<BitSet> = vec_variants
            .iter()
            .filter_map(|variant| {
                let encoded_kmer =
                    IntT::encode_kmer(&variant.sequence.get_range(0, data_info.k_graph + 1));
                kmer_2_samples.get(&encoded_kmer).cloned()
            })
            .collect();

        // compute missing samples, including `0/1` as missing
        let mut missing_samples = 0;
        let mut ref_present = false;
        let mut alt_present = false;

        for i in 0..data_info.sample_names.len() {
            let in_ref = bitset_vec[0].contains(i);
            let in_alt = bitset_vec[1].contains(i);

            if !in_ref && !in_alt {
                missing_samples += 1;
            } else if in_ref && in_alt {
                missing_samples += 1; // consider heterozygous calls as missing
            } else if in_ref {
                ref_present = true;
            } else {
                alt_present = true;
            }
        }

        let proportion_missing = missing_samples as f32 / data_info.sample_names.len() as f32;

        // filter indels based on proportion of missing data and only keep true variants
        if proportion_missing <= config.max_missing && ref_present && alt_present {
            nb_indels += 1;

            // get inserts and 1st/last k-mers
            let (vec_inserts, last_kmer) = extract_middle_bases(vec_variants, data_info.k_graph);
            let first_kmer = vec_variants[0].sequence.decode()[..data_info.k_graph].to_string();

            // determine the most frequent variant (REF) and the other (ALT)
            let mut variants: Vec<(String, usize, &BitSet)> = vec_inserts
                .iter()
                .zip(&bitset_vec)
                .map(|(seq, bitset)| (seq.clone(), bitset.len(), bitset))
                .collect();

            // sort by frequency (descending) to find the most frequent variant
            variants.sort_by(|a, b| b.1.cmp(&a.1));

            let (ref_allele, _ref_count, ref_bitset) = &variants[0]; // most frequent (REF)
            let (alt_allele, _alt_count, alt_bitset) = &variants[1]; // less frequent (ALT)

            // Generate sample genotype calls
            let sample_calls: Vec<String> = data_info
                .sample_names
                .iter()
                .enumerate()
                .map(|(i, _sample)| {
                    let in_ref = ref_bitset.contains(i);
                    let in_alt = alt_bitset.contains(i);

                    match (in_ref, in_alt) {
                        (true, true) => "0/1".to_string(), // both variants present (strain mixture)
                        (true, false) => "0".to_string(),  // only REF
                        (false, true) => "1".to_string(),  // only ALT
                        (false, false) => ".".to_string(), // missing data
                    }
                })
                .collect();

            // Write the VCF line
            writeln!(
                writer,
                ".\t.\t.\t{}\t{}\t.\tbefore={};after={}\t.\tGT\t{}",
                ref_allele,
                alt_allele,
                first_kmer,
                last_kmer,
                sample_calls.join("\t")
            )
            .unwrap();
        }
    }

    log::info!("{nb_indels} indels");

    // return entry k-mers of indels for SNP processing
    entries_indels
}

// dereplicate indel groups: choose shortest between 'forward' and 'reverse-complement';
// this is equivalent to indel realigning in read-alignment (useful in repeats)
fn dereplicate_indels<IntT: for<'a> UInt<'a>>(
    indel_groups: VariantGroups<IntT>,
    k_graph: usize,
) -> (VariantGroups<IntT>, HashSet<IntT>) {
    let mut entries_indels: HashSet<IntT> = HashSet::new();
    let mut final_indels: VariantGroups<IntT> = HashMap::new();

    // create a vector of keys and their corresponding total sequence length, and sort it in increasing order
    // we use the IntT value of the entry k-mer as tie breaker to get a stable list
    let mut sorted_extremities: Vec<((IntT, IntT), usize)> = indel_groups
        .iter()
        .map(|(key, variants)| {
            // calculate total sequence length
            let total_length: usize = variants
                .iter()
                .map(|variant| variant.sequence.decode().len())
                .sum();
            (*key, total_length)
        })
        .collect();

    sorted_extremities.sort_by(|a, b| {
        a.1.cmp(&b.1) // sort by sum of sequence lengths
            .then_with(|| a.0 .0.cmp(&b.0 .0)) // sort by the first IntT value of the key when there's a tie
    });

    for (combined_ext, _) in sorted_extremities {
        let vec_variants = indel_groups.get(&combined_ext).unwrap();
        if !entries_indels.contains(&combined_ext.0) {
            // save indel k-mers
            let rc_1 = IntT::rev_comp(combined_ext.0, k_graph);
            let rc_2 = IntT::rev_comp(combined_ext.1, k_graph);
            entries_indels.insert(combined_ext.0);
            entries_indels.insert(rc_1);
            entries_indels.insert(combined_ext.1);
            entries_indels.insert(rc_2);
            // save indel group
            final_indels.insert(combined_ext, vec_variants.clone());
        }
    }

    (final_indels, entries_indels)
}

// extract inserts of an indel group
fn extract_middle_bases(vec_variants: &[VariantInfo], k_graph: usize) -> (Vec<String>, String) {
    // collect all sequences without the first kmer
    let reduced_seq: Vec<String> = vec_variants
        .iter()
        .map(|variant| {
            let seq = variant.sequence.decode();
            seq[k_graph..].to_string()
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
    if last_kmer.len() > k_graph {
        last_kmer = last_kmer[..k_graph].to_string();
    }

    // extract 'middle-bases' (remove last kmer from reduced sequences -> only middle base left)
    let mut vec_middles: Vec<String> = Vec::new();
    for seq in &reduced_seq {
        let pos2_end = seq.len() - n_nucl;
        let mut middle_bases = &seq[..pos2_end];
        if middle_bases.is_empty() {
            middle_bases = "-";
        }
        vec_middles.push(middle_bases.to_string());
    }

    (vec_middles, last_kmer)
}
