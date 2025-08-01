//! Type for building a split k-mer dictionary from a fasta/fastq file.
//!
//! The dictionary has the concatenated left and right part of the split k-mer
//! (bit-encoded DNA bases) as keys, and ASCII IUPAC bases as values for the
//! middle base.
//!
//! Prefer to use [`crate::merge_ska_dict::build_and_merge`] over this
//! function directly.
//!
//! If you want to use this directly, build with [`SkaDict::new()`].
//!
//! These should then be converted or merged into to [`crate::merge_ska_dict::MergeSkaDict`], which
//! lets you do more useful things with them.
//!
//! Note that you should use an appropriate `sample_idx` when building if you
//! know how many you will be merging, which will let you use [`crate::merge_ska_dict::MergeSkaDict::merge()`].
//! Otherwise you will need to use the slower [`crate::merge_ska_dict::MergeSkaDict::extend()`]

use std::cmp::Ordering;

use hashbrown::HashMap;
#[cfg(not(feature = "wasm"))]
extern crate needletail;
#[cfg(not(feature = "wasm"))]
use needletail::{parse_fastx_file, parser::Format};

pub mod split_kmer;
use super::QualOpts;
use crate::ska_dict::split_kmer::SplitKmer;

pub mod bit_encoding;
use crate::ska_dict::bit_encoding::{decode_base, UInt, IUPAC};

pub mod bloom_filter;
use crate::ska_dict::bloom_filter::KmerFilter;

pub mod nthash;

#[cfg(feature = "wasm")]
use std::io::Read;
#[cfg(feature = "wasm")]
use crate::fastx_wasm::{open_fasta, open_fastq, ReaderEnum};
#[cfg(feature = "wasm")]
use seq_io::fasta::Reader as FastaReader;
#[cfg(feature = "wasm")]
use seq_io::fasta::Record as FastaRecord;
#[cfg(feature = "wasm")]
use seq_io::fastq::Reader as FastqReader;
#[cfg(feature = "wasm")]
use seq_io::fastq::Record as FastqRecord;
// #[cfg(feature = "wasm")]
// use wasm_bindgen::prelude::*;
#[cfg(feature = "wasm")]
use wasm_bindgen_file_reader::WebSysFile;

/// Holds the split-kmer dictionary, and basic information such as k-mer size.
#[derive(Debug, Clone, Default)]
pub struct SkaDict<IntT> {
    /// K-mer size
    k: usize,
    /// Whether reverse-complement was counted
    rc: bool,
    /// Sample index, if being used in a merge
    sample_idx: usize,
    /// Sample name
    name: String,
    /// Split k-mer dictionary split-k:middle-base
    split_kmers: HashMap<IntT, u8>,
    /// A bloom filter for counting from fastq files
    kmer_filter: KmerFilter,
}

impl<IntT> SkaDict<IntT>
where
    IntT: for<'a> UInt<'a>,
{
    /// Adds a split-kmer and middle base to dictionary.
    fn add_to_dict(&mut self, kmer: IntT, base: u8) {
        self.split_kmers
            .entry(kmer)
            .and_modify(|b| *b = IUPAC[base as usize * 256 + *b as usize])
            .or_insert(decode_base(base));
    }

    /// Adds a split k-mer which is a self-rc to the dict
    /// This requires amibguity of middle_base + rc(middle_base) to be added
    fn add_palindrome_to_dict(&mut self, kmer: IntT, base: u8) {
        self.split_kmers
            .entry(kmer)
            .and_modify(|b| {
                *b = match b {
                    b'W' => {
                        if base == 0 || base == 2 {
                            b'W'
                        } else {
                            b'N'
                        }
                    }
                    b'S' => {
                        if base == 0 || base == 2 {
                            b'N'
                        } else {
                            b'S'
                        }
                    }
                    b'N' => b'N',
                    _ => panic!("Palindrome middle base not W/S: {}", *b as char),
                }
            })
            .or_insert(match base {
                0 | 2 => b'W', // A or T
                1 | 3 => b'S', // C or G
                _ => panic!("Base encoding error: {}", base as char),
            });
    }

    /// Iterates through all the k-mers from an input fastx file and adds them
    /// to the dictionary
    #[cfg(not(feature = "wasm"))]
    fn add_file_kmers(
        &mut self,
        filename: &str,
        is_reads: bool,
        qual: &QualOpts,
        proportion_reads: Option<f64>,
    ) {
        let mut step: usize = 1;

        if let Some(prop) = proportion_reads {
            step = (1.0 / prop).round() as usize;
        }

        let mut reader =
            parse_fastx_file(filename).unwrap_or_else(|_| panic!("Invalid path/file: {filename}"));

        let mut iter_reads: usize = 0;
        while let Some(record) = reader.next() {
            if !iter_reads.is_multiple_of(step) {
                iter_reads += 1;
                continue;
            } else {
                iter_reads += 1;
            }

            let seqrec = record.expect("Invalid FASTA/Q record");
            let kmer_opt = SplitKmer::new(
                seqrec.seq(),
                seqrec.num_bases(),
                seqrec.qual(),
                self.k,
                self.rc,
                qual.min_qual,
                qual.qual_filter,
                is_reads,
            );
            if let Some(mut kmer_it) = kmer_opt {
                if !is_reads
                    || (kmer_it.middle_base_qual()
                        && Ordering::is_eq(self.kmer_filter.filter(&kmer_it)))
                {
                    let (kmer, base, _rc) = kmer_it.get_curr_kmer();
                    if kmer_it.self_palindrome() {
                        self.add_palindrome_to_dict(kmer, base);
                    } else {
                        self.add_to_dict(kmer, base);
                    }
                }
                while let Some((kmer, base, _rc)) = kmer_it.get_next_kmer() {
                    if !is_reads
                        || (kmer_it.middle_base_qual()
                            && Ordering::is_eq(self.kmer_filter.filter(&kmer_it)))
                    {
                        if kmer_it.self_palindrome() {
                            self.add_palindrome_to_dict(kmer, base);
                        } else {
                            self.add_to_dict(kmer, base);
                        }
                    }
                }
            }
        }
    }

    /// Iterates through all the k-mers from an input fastx file and adds them
    /// to the dictionary
    #[cfg(feature = "wasm")]
    pub fn add_file_kmers<F: Read>(
        &mut self,
        file: &mut F,
        is_reads: bool,
        qual: &QualOpts,
        proportion_reads: Option<f64>,
    ) {

        enum ReaderType<'a, F: Read + 'a> {
            Fasta(FastaReader<ReaderEnum<'a, F>>), // replace with actual type
            Fastq(FastqReader<ReaderEnum<'a, F>>), // replace with actual type
        }

        let mut reader;

        if !is_reads {
            reader = ReaderType::Fasta(open_fasta(file));
        } else {
            reader = ReaderType::Fastq(open_fastq(file));
        }

        let step = 1 as f64 / proportion_reads.unwrap();

        let mut iter_reads = 0;
        while let Some((seq, seq_len)) = match reader {
            ReaderType::Fasta(ref mut r) => {
                if let Some(record) = r.next(){
                    let seqrec = record.expect("Invalid FASTA record");
                    // There can be \n in the sequence, its ascii code is 10
                    let seq: Vec<u8> = seqrec.seq().to_vec().iter().filter(|&x| *x != 10).cloned().collect();
                    let seq_len = seq.len();
                    Some((seq, seq_len))
                } else {
                    None
                }
            }
            ReaderType::Fastq(ref mut r) => {
                if let Some(record) = r.next(){
                    let seqrec = record.expect("Invalid FASTQ record");
                    // There can be \n in the sequence, its ascii code is 10
                    let seq: Vec<u8> = seqrec.seq().to_vec().iter().filter(|&x| *x != 10).cloned().collect();
                    let seq_len = seq.len();
                    Some((seq, seq_len))
                } else {
                    None
                }
            }
        } {
            if (iter_reads as f64 % step) as i32 != 0 {
                iter_reads += 1;
                continue;
            } else {
                iter_reads += 1;
            }

            let kmer_opt = SplitKmer::new(
                seq.into(),
                seq_len,
                None,
                self.k,
                self.rc,
                qual.min_qual,
                qual.qual_filter,
                is_reads,
            );
            if let Some(mut kmer_it) = kmer_opt {
                if !is_reads
                    || (kmer_it.middle_base_qual()
                        && Ordering::is_eq(self.kmer_filter.filter(&kmer_it)))
                {
                    let (kmer, base, _rc) = kmer_it.get_curr_kmer();
                    if kmer_it.self_palindrome() {
                        self.add_palindrome_to_dict(kmer, base);
                    } else {
                        self.add_to_dict(kmer, base);
                    }
                }
                while let Some((kmer, base, _rc)) = kmer_it.get_next_kmer() {
                    if !is_reads
                        || (kmer_it.middle_base_qual()
                            && Ordering::is_eq(self.kmer_filter.filter(&kmer_it)))
                    {
                        if kmer_it.self_palindrome() {
                            self.add_palindrome_to_dict(kmer, base);
                        } else {
                            self.add_to_dict(kmer, base);
                        }
                    }
                }
            }
        }
    }

    /// Build a split-kmer dictionary from input fastx file(s)
    ///
    /// Prefer to use [`crate::merge_ska_dict::build_and_merge()`] over this
    /// function directly.
    ///
    /// # Examples
    ///
    /// To build with a FASTA
    /// ```
    /// use ska::ska_dict::SkaDict;
    /// use ska::{QualOpts, QualFilter};
    ///
    /// let k = 31;
    /// let sample_idx = 0;
    /// let quality = QualOpts {min_count: 1, min_qual: 0, qual_filter: QualFilter::NoFilter};
    /// let proportion_reads = None;
    /// let ska_dict = SkaDict::<u64>::new(k, sample_idx, (&"tests/test_files_in/test_1.fa", None), "test_1", true, &quality, proportion_reads);
    /// ```
    ///
    /// With FASTQ pair, only allowing k-mers with a count over 2, and where all
    /// bases have a PHRED score of 20 or more
    /// ```
    /// use ska::ska_dict::SkaDict;
    /// use ska::{QualOpts, QualFilter};
    ///
    /// let quality = QualOpts {min_count: 2, min_qual: 20, qual_filter: QualFilter::Middle};
    /// let k = 9;
    /// let sample_idx = 0;
    /// let proportion_reads = None;
    /// let ska_dict = SkaDict::<u64>::new(k, sample_idx,
    ///                             (&"tests/test_files_in/test_1_fwd.fastq.gz",
    ///                             Some(&"tests/test_files_in/test_2_fwd.fastq.gz".to_string())),
    ///                             "sample1", true, &quality, proportion_reads);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if:
    /// - k-mer length is invalid (<5, >63 or even)
    /// - Input file cannot be read
    /// - Input file contains invalid fastx record
    /// - Input file contains no valid sequence to find at least on split k-mer
    #[cfg(not(feature = "wasm"))]
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        k: usize,
        sample_idx: usize,
        files: (&str, Option<&String>),
        name: &str,
        rc: bool,
        qual: &QualOpts,
        proportion_reads: Option<f64>,
    ) -> Self {
        if !(5..=63).contains(&k) || k.is_multiple_of(2) {
            panic!("Invalid k-mer length");
        }

        let mut sk_dict = Self {
            k,
            rc,
            sample_idx,
            name: name.to_string(),
            split_kmers: HashMap::default(),
            kmer_filter: KmerFilter::new(qual.min_count),
        };

        // Check if we're working with reads, and initalise the CM filter if so
        let mut reader_peek =
            parse_fastx_file(files.0).unwrap_or_else(|_| panic!("Invalid path/file: {}", files.0));
        let seq_peek = reader_peek
            .next()
            .expect("Invalid FASTA/Q record")
            .expect("Invalid FASTA/Q record");
        let mut is_reads = false;
        if seq_peek.format() == Format::Fastq {
            sk_dict.kmer_filter.init();
            is_reads = true;
        }

        // Build the dict
        sk_dict.add_file_kmers(files.0, is_reads, qual, proportion_reads);
        if let Some(second_filename) = files.1 {
            sk_dict.add_file_kmers(second_filename, is_reads, qual, proportion_reads);
        }

        if sk_dict.ksize() == 0 {
            panic!("{} has no valid sequence", files.0);
        }
        sk_dict
    }

    /// Build a split-kmer dictionary from input fastx file(s)
    ///
    /// Prefer to use [`crate::merge_ska_dict::build_and_merge()`] over this
    /// function directly.
    ///
    /// # Examples
    ///
    /// To build with a FASTA
    /// ```
    /// use ska::ska_dict::SkaDict;
    /// use ska::{QualOpts, QualFilter};
    ///
    /// let k = 31;
    /// let sample_idx = 0;
    /// let quality = QualOpts {min_count: 1, min_qual: 0, qual_filter: QualFilter::NoFilter};
    /// let ska_dict = SkaDict::<u64>::new(k, sample_idx, (&"tests/test_files_in/test_1.fa", None), "test_1", true, &quality);
    /// ```
    ///
    /// With FASTQ pair, only allowing k-mers with a count over 2, and where all
    /// bases have a PHRED score of 20 or more
    /// ```
    /// use ska::ska_dict::SkaDict;
    /// use ska::{QualOpts, QualFilter};
    ///
    /// let quality = QualOpts {min_count: 2, min_qual: 20, qual_filter: QualFilter::Middle};
    /// let k = 9;
    /// let sample_idx = 0;
    /// let ska_dict = SkaDict::<u64>::new(k, sample_idx,
    ///                             (&"tests/test_files_in/test_1_fwd.fastq.gz",
    ///                             Some(&"tests/test_files_in/test_2_fwd.fastq.gz".to_string())),
    ///                             "sample1", true, &quality);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if:
    /// - k-mer length is invalid (<5, >63 or even)
    /// - Input file cannot be read
    /// - Input file contains invalid fastx record
    /// - Input file contains no valid sequence to find at least on split k-mer
    #[cfg(feature = "wasm")]
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        k: usize,
        sample_idx: usize,
        input_files: (&web_sys::File, Option<&web_sys::File>),
        name: &str,
        rc: bool,
        qual: &QualOpts,
        proportion_reads: Option<f64>,
    ) -> Self {
        if !(5..=63).contains(&k) || k % 2 == 0 {
            panic!("Invalid k-mer length");
        }

        let file_name = input_files.0.name();
        let mut file_type = file_name.split('.').nth(file_name.split('.').count() - 1).unwrap();
        if file_type == "gz" {
            file_type = file_name.split('.').nth(file_name.split('.').count() - 2).unwrap();
        }

        let is_reads : bool;
        if ["fasta", "fa"].contains(&file_type) {
            is_reads = false;
        } else if ["fastq", "fq"].contains(&file_type) {
            is_reads = true;
        } else {
            panic!("Unsupported file type.")
        }

        let mut sk_dict = Self {
            k,
            rc,
            sample_idx,
            name: name.to_string(),
            split_kmers: HashMap::default(),
            kmer_filter: KmerFilter::new(qual.min_count),
        };

        // TODO Check if we're working with reads, and initalise the CM filter if so
        // and self.is_reads
        /*
        let mut reader_peek =
            parse_fastx_file(files.0).unwrap_or_else(|_| panic!("Invalid path/file: {}", files.0));
        let seq_peek = reader_peek
            .next()
            .expect("Invalid FASTA/Q record")
            .expect("Invalid FASTA/Q record");
        let mut is_reads = false;
        if seq_peek.format() == Format::Fastq {
            sk_dict.kmer_filter.init();
            is_reads = true;
        }
        */

        // Build the dict
        sk_dict.add_file_kmers(&mut WebSysFile::new(input_files.0.clone()), is_reads, qual, proportion_reads);

        if let Some(second_file) = input_files.1 {
            sk_dict.add_file_kmers(&mut WebSysFile::new(second_file.clone()), is_reads, qual, proportion_reads);
        }

        if sk_dict.ksize() == 0 {
            panic!("File has no valid sequence");
        }
        sk_dict
    }

    /// K-mer length used for split k-mers
    pub fn kmer_len(&self) -> usize {
        self.k
    }

    /// Whether reverse-complement was counted
    pub fn rc(&self) -> bool {
        self.rc
    }

    /// Total number of split k-mers in the dictionary
    pub fn ksize(&self) -> usize {
        self.split_kmers.len()
    }

    /// Sample index for merging
    pub fn idx(&self) -> usize {
        self.sample_idx
    }

    /// Split k-mer dictonary
    pub fn kmers(&self) -> &HashMap<IntT, u8> {
        &self.split_kmers
    }

    /// Sample name
    pub fn name(&self) -> &String {
        &self.name
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_add_palindrome_to_dict() {
        // Initialize the test object
        let mut test_obj = SkaDict::<u64>::default();

        // Test case 1: Updating existing entry
        test_obj.split_kmers.insert(123, b'W');
        test_obj.add_palindrome_to_dict(123, 1);
        assert_eq!(test_obj.split_kmers[&123], b'N');

        // Test case 2: Adding new entry with base 0
        test_obj.split_kmers.clear();
        test_obj.add_palindrome_to_dict(456, 0);
        assert_eq!(test_obj.split_kmers[&456], b'W');

        // Test case 3: Adding new entry with base 3
        test_obj.split_kmers.clear();
        test_obj.add_palindrome_to_dict(789, 3);
        assert_eq!(test_obj.split_kmers[&789], b'S');

        // Test case 4: Updating existing twice
        test_obj.split_kmers.insert(123, b'W');
        test_obj.add_palindrome_to_dict(123, 1);
        test_obj.add_palindrome_to_dict(123, 1);
        assert_eq!(test_obj.split_kmers[&123], b'N');
    }

    #[test]
    #[should_panic]
    fn test_panic_add_palindrome_to_dict() {
        // Test case 4: Panicking with invalid base
        let mut test_obj_panic = SkaDict::<u64>::default();
        test_obj_panic.add_palindrome_to_dict(987, 5);
    }

    #[test]
    #[should_panic]
    fn test_panic2_add_palindrome_to_dict() {
        // Test case 5: Panicking with invalid middle base
        let mut test_obj_panic = SkaDict::<u64>::default();
        test_obj_panic.split_kmers.clear();
        test_obj_panic.split_kmers.insert(555, b'A');
        test_obj_panic.add_palindrome_to_dict(555, 1);
    }
}
