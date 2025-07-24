//! Functions for writing output in ska lo
use hashbrown::{HashMap, HashSet};
use std::fs::File;
use std::io::Write;

use crate::skalo::utils::Config;

/// the function outputs a SNP alignment, and a VCF file and a pseudo-genome alignment
/// if a reference genome has been provided
pub fn create_fasta_and_vcf(
    genome_name: String,
    mut genome_seq: Vec<u8>,
    sample_names: Vec<String>,
    variant_map: HashMap<u32, Vec<char>>,
    config: &Config,
) {
    // replace non-ATGCN characters with 'N' in genome_seq
    for base in genome_seq.iter_mut() {
        match *base as char {
            'A' | 'T' | 'G' | 'C' | 'N' => {} // Valid bases remain unchanged
            _ => *base = b'N',                // Replace other characters with 'N'
        }
    }

    // sort variants by positions (increasing u32 key)
    let mut sorted_map: Vec<_> = variant_map.into_iter().collect();
    sorted_map.sort_by_key(|&(key, _)| key);

    let mut sequences: Vec<String> = vec![String::new(); sample_names.len()];
    let mut genome_alignments: Option<Vec<String>> = None; // only create if genome_seq is not empty
    let mut vcf_records: Vec<(u32, char, Vec<char>)> = Vec::new();

    if !genome_seq.is_empty() {
        genome_alignments = Some(vec![String::new(); sample_names.len()]);
    }

    let mut current_snp_index = 0;

    // iterate through genome positions if genome_seq is provided
    let genome_length = if !genome_seq.is_empty() {
        genome_seq.len() as u32
    } else {
        sorted_map.last().map_or(0, |(pos, _)| *pos + 1)
    };

    for pos in 0..genome_length {
        if current_snp_index < sorted_map.len() && sorted_map[current_snp_index].0 == pos {
            // SNP position
            let (snp_pos, vec_chars) = &sorted_map[current_snp_index];

            if let Some(ref_genome_alignments) = genome_alignments.as_mut() {
                let reference_base = genome_seq[*snp_pos as usize] as char;

                // push the variant data into vcf_records
                vcf_records.push((*snp_pos, reference_base, vec_chars.clone()));

                // build genome alignments
                for (i, &char) in vec_chars.iter().enumerate() {
                    ref_genome_alignments[i].push(char);
                }
            }

            // build SNP-only sequences
            for (i, &char) in vec_chars.iter().enumerate() {
                sequences[i].push(char);
            }

            current_snp_index += 1;
        } else if let Some(ref_genome_alignments) = genome_alignments.as_mut() {
            // non-SNP position: add reference base to genome alignment
            let ref_base = genome_seq[pos as usize] as char;
            for genome_alignment in ref_genome_alignments.iter_mut() {
                genome_alignment.push(ref_base);
            }
        }
    }

    // write SNP alignment in FASTA format
    let snp_filename = format!("{}_snps.fas", config.output_name);
    let mut snp_output = File::create(snp_filename).expect("Unable to create SNP file");
    for (name, sequence) in sample_names.iter().zip(sequences.iter()) {
        writeln!(snp_output, ">{name}").expect("Error writing to SNP file");
        writeln!(snp_output, "{sequence}").expect("Error writing to SNP file");
    }

    if !genome_seq.is_empty() {
        // write pseudo-genomes in FASTA format
        let genome_filename = format!("{}_pseudo_genomes.fas", config.output_name);
        let genome_output =
            File::create(genome_filename).expect("Unable to create genome alignment file");
        if let Some(ref_genome_alignments) = genome_alignments {
            for (name, alignment) in sample_names.iter().zip(ref_genome_alignments.iter()) {
                writeln!(&genome_output, ">{name}")
                    .expect("Error writing to genome alignment file");
                writeln!(&genome_output, "{alignment}")
                    .expect("Error writing to genome alignment file");
            }
        }

        // write variants in VCF format
        let vcf_filename = format!("{}_snps.vcf", config.output_name);
        let mut vcf_output = File::create(vcf_filename).expect("Unable to create VCF file");
        writeln!(vcf_output, "##fileformat=VCFv4.2").expect("Error writing VCF header");
        writeln!(
            vcf_output,
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}",
            sample_names.join("\t")
        )
        .expect("Error writing VCF header");

        for (pos, reference_base, vec_chars) in vcf_records {
            let alt_bases: Vec<char> = vec_chars
                .iter()
                .cloned()
                .filter(|&c| c != reference_base && c != '-' && c != 'N')
                .collect::<HashSet<_>>() // deduplicate alternative bases
                .into_iter()
                .collect();

            let genotypes: Vec<String> = vec_chars
                .iter()
                .map(|&c| {
                    if c == reference_base {
                        "0".to_string()
                    } else if c == '-' || c == 'N' {
                        ".".to_string() // missing or ambiguous data
                    } else if let Some(alt_index) = alt_bases.iter().position(|&alt| alt == c) {
                        (alt_index + 1).to_string() // ALT indices in VCF are 1-based
                    } else {
                        ".".to_string() // unexpected case
                    }
                })
                .collect();

            writeln!(
                vcf_output,
                "{}\t{}\t.\t{}\t{}\t.\t.\t.\tGT\t{}",
                genome_name,
                pos + 1, // VCF positions are 1-based
                reference_base,
                alt_bases
                    .iter()
                    .map(|&c| c.to_string())
                    .collect::<Vec<_>>()
                    .join(","),
                genotypes.join("\t")
            )
            .expect("Error writing VCF record");
        }
    }
}
