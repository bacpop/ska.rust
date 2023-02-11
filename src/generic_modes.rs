use crate::cli::{FileType, FilterType};
use crate::io_utils::set_ostream;
use crate::merge_ska_dict::MergeSkaDict;
use crate::merge_ska_array::MergeSkaArray;
use crate::ska_dict::bit_encoding::RevComp;
use crate::ska_ref::RefSka;

pub fn align<IntT: for<'a> RevComp<'a>>(
    ska_array: &mut MergeSkaArray<IntT>,
    output: &Option<String>,
    filter: &FilterType,
    min_freq: f64,
) {
    // In debug mode (cannot be set from CLI, give details)
    log::debug!("{ska_array}");

    // Apply filters
    let update_kmers = false;
    let filter_threshold = f64::ceil(ska_array.nsamples() as f64 * min_freq) as usize;
    log::info!("Applying filters: threshold={filter_threshold} constant_site_filter={filter}");
    ska_array.filter(filter_threshold, filter, update_kmers);

    // Write out to file/stdout
    log::info!("Writing alignment");
    let mut out_stream = set_ostream(output);
    ska_array
        .write_fasta(&mut out_stream)
        .expect("Couldn't write output fasta");
}

pub fn map<IntT: for<'a> RevComp<'a>>(
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

pub fn merge<IntT: for<'a> RevComp<'a>>(
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

pub fn delete<IntT: for<'a> RevComp<'a>>(
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

pub fn weed<IntT: for<'a> RevComp<'a>>(
    ska_array: &mut MergeSkaArray<IntT>,
    weed_file: &Option<String>,
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

        log::info!("Removing weed k-mers");
        ska_array.weed(&ska_weed);
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

pub fn save_skf<IntT: for<'a> RevComp<'a>>(ska_dict: &MergeSkaDict<IntT>, out_file: &str) {
    log::info!("Converting to array representation and saving");
    let ska_array = MergeSkaArray::new(ska_dict);
    ska_array
        .save(out_file)
        .expect("Failed to save output file");
}
