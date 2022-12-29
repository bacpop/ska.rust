// Class for split-kmers from one input file

extern crate needletail;
use hashbrown::HashMap;
use needletail::{parse_fastx_file, parser::Format};

pub mod split_kmer;
use crate::ska_dict::split_kmer::SplitKmer;

pub mod bit_encoding;
use crate::ska_dict::bit_encoding::{decode_base, IUPAC};

pub mod count_min_filter;
use crate::ska_dict::count_min_filter::CountMin;
const CM_WIDTH: usize = 134217728; // 2^27
const CM_HEIGHT: usize = 4;

pub struct SkaDict {
    k: usize,
    rc: bool,
    sample_idx: usize,
    name: String,
    split_kmers: HashMap<u64, u8>,
    cm_filter: CountMin,
}

impl SkaDict {
    fn add_to_dict(&mut self, kmer: u64, base: u8, is_reads: bool) {
        if !is_reads || self.cm_filter.filter(kmer) {
            self.split_kmers
                .entry(kmer)
                .and_modify(|b| *b = IUPAC[base as usize * 256 + *b as usize])
                .or_insert(decode_base(base));
        }
    }

    // Iterates through all the k-mers from a file
    fn add_file_kmers(&mut self, filename: &str, is_reads: bool, min_qual: u8) {
        let mut reader =
            parse_fastx_file(filename).expect(&format!("Invalid path/file: {}", filename));
        while let Some(record) = reader.next() {
            let seqrec = record.expect("Invalid FASTA/Q record");
            let kmer_opt = SplitKmer::new(seqrec.seq(), seqrec.num_bases(), seqrec.qual(), self.k, self.rc, min_qual);
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

    pub fn kmer_len(&self) -> usize {
        self.k
    }

    pub fn rc(&self) -> bool {
        self.rc
    }

    pub fn ksize(&self) -> usize {
        self.split_kmers.len()
    }

    pub fn idx(&self) -> usize {
        self.sample_idx
    }

    pub fn kmers(&self) -> &HashMap<u64, u8> {
        &self.split_kmers
    }

    pub fn name(&self) -> &String {
        &self.name
    }

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
}
