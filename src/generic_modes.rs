//! Main control of most CLI functions, generic over `u64` and `u128`.
//!
//! Includes functions for `ska align`, `ska map`, `ska merge`, `ska distance`, `ska delete`,
//! and `ska weed`. These are needed as when loading an skf file we don't know
//! what the int type used is, so we want to then dispatch the sucessfully loaded
//! file to a generic function.

use std::io::Write;

use crate::cli::{FileType, FilterType};
use crate::io_utils::set_ostream;
use crate::merge_ska_array::MergeSkaArray;
use crate::merge_ska_dict::MergeSkaDict;
use crate::ska_dict::bit_encoding::UInt;
use crate::ska_ref::RefSka;

/// Filters alignment, and prints it out
pub fn align<IntT: for<'a> UInt<'a>>(
    ska_array: &mut MergeSkaArray<IntT>,
    output: &Option<String>,
    filter: &FilterType,
    min_freq: f64,
) {
    // In debug mode (cannot be set from CLI, give details)
    log::debug!("{ska_array}");

    // Apply filters
    apply_filters(ska_array, min_freq, filter);

    // Write out to file/stdout
    log::info!("Writing alignment");
    let mut out_stream = set_ostream(output);
    ska_array
        .write_fasta(&mut out_stream)
        .expect("Couldn't write output fasta");
}

/// Code for `ska map`
///
/// Convert array to dictionary representation, runs map against a reference,
/// prints out in requested format.
pub fn map<IntT: for<'a> UInt<'a>>(
    ska_array: &mut MergeSkaArray<IntT>,
    ska_ref: &mut RefSka<IntT>,
    output: &Option<String>,
    format: &FileType,
    threads: usize,
) {
    log::info!("Converting skf to dictionary");
    let ska_dict = ska_array.to_dict();

    log::info!("Mapping");
    ska_ref.map(&ska_dict);

    let mut out_stream = set_ostream(output);
    match format {
        FileType::Aln => {
            log::info!("Generating alignment with {} threads", threads);
            ska_ref
                .write_aln(&mut out_stream, threads)
                .expect("Failed to write output alignment");
        }
        FileType::Vcf => {
            log::info!("Generating VCF with {} threads", threads);
            ska_ref
                .write_vcf(&mut out_stream, threads)
                .expect("Failed to write output VCF");
        }
    }
}

/// Merge multiple skf files. Need to give the first file separately to define
/// the integer type.
///
/// Subsequent files are most easily passed as a slice with `[1..]`
pub fn merge<IntT: for<'a> UInt<'a>>(
    first_array: &mut MergeSkaArray<IntT>,
    skf_files: &[String],
    output: &str,
) {
    let mut merged_dict = first_array.to_dict();

    for (file_idx, skf_in) in skf_files.iter().enumerate() {
        log::info!("Merging alignment {}", format!("{}", file_idx + 1));
        let next_array = MergeSkaArray::load(skf_in)
            .expect("Failed to load input file (inconsistent k-mer lengths?)");
        merged_dict.extend(&mut next_array.to_dict());
    }

    log::info!("Converting and saving merged alignment");
    let merged_array = MergeSkaArray::new(&merged_dict);
    merged_array
        .save(format!("{output}.skf").as_str())
        .expect("Failed to save output file");
}

/// Apply a constant/ambigbuous filter [`FilterType`] and a minimum frequency
/// filter to the array of variants, updating the object.
///
/// Returns the number of removed sites
pub fn apply_filters<IntT: for<'a> UInt<'a>>(
    ska_array: &mut MergeSkaArray<IntT>,
    min_freq: f64,
    filter: &FilterType,
) -> i32 {
    let update_kmers = false;
    let filter_threshold = f64::ceil(ska_array.nsamples() as f64 * min_freq) as usize;
    log::info!("Applying filters: threshold={filter_threshold} constant_site_filter={filter}");
    ska_array.filter(filter_threshold, filter, update_kmers)
}

/// Calculate distances between samples
///
/// Also applies filters and prints them out in 'long' form
pub fn distance<IntT: for<'a> UInt<'a>>(
    ska_array: &mut MergeSkaArray<IntT>,
    output_prefix: &Option<String>,
    min_freq: f64,
    filt_ambig: bool,
    threads: usize,
) {
    // In debug mode (cannot be set from CLI, give details)
    log::debug!("{ska_array}");

    let constant = apply_filters(ska_array, min_freq, &FilterType::NoConst);
    if filt_ambig || (min_freq * ska_array.nsamples() as f64 >= 1.0) {
        apply_filters(ska_array, min_freq, &FilterType::NoAmbigOrConst);
    }

    log::info!("Calculating distances");
    if threads > 1 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .unwrap();
    }
    let distances = ska_array.distance(constant as f64);

    // Write out the distances (long form)
    let mut f = set_ostream(output_prefix);
    writeln!(&mut f, "Sample1\tSample2\tDistance\tMismatches").unwrap();
    let sample_names = ska_array.names();
    for (idx, (dist_vec, sample1)) in distances.iter().zip(sample_names).enumerate() {
        for ((distance, mismatches), j) in dist_vec.iter().zip(std::ops::Range {
            start: (idx + 1),
            end: sample_names.len(),
        }) {
            writeln!(
                &mut f,
                "{sample1}\t{}\t{distance:.2}\t{mismatches:.5}",
                sample_names[j]
            )
            .unwrap();
        }
    }
}

/// Delete files with passed names in the given array
pub fn delete<IntT: for<'a> UInt<'a>>(
    ska_array: &mut MergeSkaArray<IntT>,
    delete_names: &[&str],
    out_file: &str,
) {
    log::info!("Deleting samples");
    ska_array.delete_samples(delete_names);

    log::info!("Saving modified skf file");
    ska_array
        .save(out_file)
        .expect("Could not save modified array");
}

/// Remove k-mers, and optionally apply filters to an array
pub fn weed<IntT: for<'a> UInt<'a>>(
    ska_array: &mut MergeSkaArray<IntT>,
    weed_file: &Option<String>,
    reverse: bool,
    min_freq: f64,
    filter: &FilterType,
    out_file: &str,
) {
    if let Some(weed_fasta) = weed_file {
        log::info!(
            "Making skf of weed file k={} rc={}",
            ska_array.kmer_len(),
            ska_array.rc()
        );
        let ska_weed = RefSka::new(ska_array.kmer_len(), weed_fasta, ska_array.rc());

        if !reverse {
            log::info!("Removing weed k-mers");
        } else {
            log::info!("Keeping only weed k-mers");
        }
        ska_array.weed(&ska_weed, reverse);
    }

    let filter_threshold = f64::floor(ska_array.nsamples() as f64 * min_freq) as usize;
    if filter_threshold > 0 || *filter != FilterType::NoFilter {
        log::info!("Applying filters: threshold={filter_threshold} constant_site_filter={filter}");
        let update_kmers = true;
        ska_array.filter(filter_threshold, filter, update_kmers);
    }

    log::info!("Saving modified skf file");
    ska_array
        .save(out_file)
        .expect("Failed to save output file");
}

/// Covnvert a dictionary representation to an array and save to file
pub fn save_skf<IntT: for<'a> UInt<'a>>(ska_dict: &MergeSkaDict<IntT>, out_file: &str) {
    log::info!("Converting to array representation and saving");
    let ska_array = MergeSkaArray::new(ska_dict);
    ska_array
        .save(out_file)
        .expect("Failed to save output file");
}
