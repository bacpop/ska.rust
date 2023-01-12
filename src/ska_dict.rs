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

extern crate needletail;
use hashbrown::HashMap;
use needletail::{parse_fastx_file, parser::Format};

pub mod split_kmer;
use crate::ska_dict::split_kmer::SplitKmer;

pub mod bit_encoding;
use crate::ska_dict::bit_encoding::{decode_base, IUPAC};

pub mod count_min_filter;
use crate::ska_dict::count_min_filter::CountMin;
const CM_WIDTH: usize = 1 << 27; // 2^27 =~ 130M
const CM_HEIGHT: usize = 4;

/// Holds the split-kmer dictionary, and basic information such as k-mer size.
#[derive(Debug, Clone)]
pub struct SkaDict {
    /// K-mer size
    k: usize,
    /// Whether reverse-complement was counted
    rc: bool,
    /// Sample index, if being used in a merge
    sample_idx: usize,
    /// Sample name
    name: String,
    /// Split k-mer dictionary split-k:middle-base
    split_kmers: HashMap<u64, u8>,
    /// A countmin filter for counting from fastq files
    cm_filter: CountMin,
}

impl SkaDict {
    /// Adds a split-kmer and middle base to dictionary. If `is_reads` then
    /// only adds if passing through the countmin filter
    fn add_to_dict(&mut self, kmer: u64, base: u8, is_reads: bool) {
        if !is_reads || self.cm_filter.filter(kmer) {
            self.split_kmers
                .entry(kmer)
                .and_modify(|b| *b = IUPAC[base as usize * 256 + *b as usize])
                .or_insert(decode_base(base));
        }
    }

    /// Iterates through all the k-mers from an input fastx file and adds them
    /// to the dictionary
    fn add_file_kmers(&mut self, filename: &str, is_reads: bool, min_qual: u8) {
        let mut reader =
            parse_fastx_file(filename).expect(&format!("Invalid path/file: {}", filename));
        while let Some(record) = reader.next() {
            let seqrec = record.expect("Invalid FASTA/Q record");
            let kmer_opt = SplitKmer::new(
                seqrec.seq(),
                seqrec.num_bases(),
                seqrec.qual(),
                self.k,
                self.rc,
                min_qual,
            );
            if kmer_opt.is_some() {
                let mut kmer_it = kmer_opt.unwrap();
                let (kmer, base, _rc) = kmer_it.get_curr_kmer();
                self.add_to_dict(kmer, base, is_reads);
                while let Some((kmer, base, _rc)) = kmer_it.get_next_kmer() {
                    self.add_to_dict(kmer, base, is_reads);
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
    /// ```rust
    /// let k = 31;
    /// let sample_idx = 0;
    /// let ska_dict = SkaDict::new(k, sample_idx, &"file.fa", None, "sample1", true, 0, 0);
    /// ```
    ///
    /// With FASTQ pair, only allowing k-mers with a count over 20, and where all
    /// bases have a PHRED score of 30 or more
    /// ```
    /// let ska_dict = SkaDict::new(k, sample_idx, &"reads_1.fastq.gz", &Some("reads_2.fastq.gz"), "sample1", true, 20, 30);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if:
    /// - k-mer length is invalid (<5, >31 or even)
    /// - Input file cannot be read
    /// - Input file contains invalid fastx record
    /// - Input file contains no valid sequence to find at least on split k-mer
    pub fn new(
        k: usize,
        sample_idx: usize,
        filename: &str,
        rev_file: &Option<String>,
        name: &str,
        rc: bool,
        min_count: u16,
        min_qual: u8,
    ) -> Self {
        if k < 5 || k > 31 || k % 2 == 0 {
            panic!("Invalid k-mer length");
        }
        // Default/empty structs
        let name = name.to_string();
        let split_kmers: HashMap<u64, u8> = HashMap::default();
        let cm_filter = CountMin::empty(CM_WIDTH, CM_HEIGHT, min_count);
        let mut sk_dict = Self {
            k,
            rc,
            sample_idx,
            name,
            split_kmers,
            cm_filter,
        };

        // Check if we're working with reads, and initalise the CM filter if so
        let mut reader_peek =
            parse_fastx_file(&filename).expect(&format!("Invalid path/file: {}", filename));
        let seq_peek = reader_peek
            .next()
            .expect("Invalid FASTA/Q record")
            .expect("Invalid FASTA/Q record");
        let mut is_reads = false;
        if seq_peek.format() == Format::Fastq {
            sk_dict.cm_filter.init();
            is_reads = true;
        }

        // Build the dict
        sk_dict.add_file_kmers(filename, is_reads, min_qual);
        if let Some(second_filename) = rev_file {
            sk_dict.add_file_kmers(second_filename, is_reads, min_qual);
        }

        if sk_dict.ksize() == 0 {
            panic!("{} has no valid sequence", filename);
        }
        return sk_dict;
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
    pub fn kmers(&self) -> &HashMap<u64, u8> {
        &self.split_kmers
    }

    /// Sample name
    pub fn name(&self) -> &String {
        &self.name
    }
}
