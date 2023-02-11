use crate::cli::FilterType;
use crate::io_utils::set_ostream;
use crate::merge_ska_array::MergeSkaArray;
use crate::ska_dict::bit_encoding::RevComp;

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
