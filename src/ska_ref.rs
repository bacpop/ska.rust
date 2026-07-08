//! Types for listing split-kmers from one input file (usually a FASTA reference sequence).
//!
//! In [`RefSka`] split k-mers are stored in a [`Vec`] so fast iteration but not fast lookup
//! is supported.
//! This generally will be a reference file, but may also be k-mers to be
//! tested/removed.
//!
//! Primarily used to support `ska map`. Functions to write out a mapped alignment
//! (multi-FASTA using [`needletail`]) or VCF (using [`noodles_vcf`]) are included.
//!
//! # Examples
//!  ```
//! use ska::ska_ref::RefSka;
//! use ska::io_utils::{load_array, set_ostream};
//!
//! // Load a saved array from build, convert to dict representation
//! let threads = 1;
//! // If you don't know whether u64 or u128, try one and handle the error
//! // see lib.rs for examples
//! let ska_dict = load_array::<u64>(&["tests/test_files_in/merge.skf".to_string()], threads).unwrap().to_dict();
//!
//! // Index a reference sequence
//! let mask_repeats = false;
//! let mask_ambiguous = false;
//! let mut ref_kmers = RefSka::new(ska_dict.kmer_len(), &"tests/test_files_in/test_ref.fa", ska_dict.rc(), mask_repeats, mask_ambiguous);
//!
//! // Run mapping, output an alignment to stdout
//! ref_kmers.map(&ska_dict);
//! let mut out_stream = set_ostream(&None);
//! ref_kmers.write_aln(&mut out_stream, threads);
//! ```

#[cfg(not(target_family = "wasm"))]
use std::io::Write;
#[cfg(not(target_family = "wasm"))]
use std::str;

use hashbrown::hash_set::Entry::*;
use hashbrown::HashSet;
#[cfg(target_family = "wasm")]
use hashbrown::HashMap;
use rayon::prelude::*;

#[cfg(not(target_family = "wasm"))]
use noodles_vcf::{
    self as vcf,
    header::record::value::{map::Contig, Map},
    record::{
        alternate_bases::Allele,
        genotypes::{keys::key, sample::Value, Keys},
        reference_bases::Base,
        AlternateBases, Genotypes, Position,
    },
};

#[cfg(not(target_family = "wasm"))]
extern crate needletail;
#[cfg(not(target_family = "wasm"))]
use ndarray::{s, Array2, ArrayView};
#[cfg(not(target_family = "wasm"))]
use needletail::{
    parse_fastx_file,
    parser::{write_fasta, Format},
};

use super::QualFilter;
pub mod aln_writer;
use crate::ska_ref::aln_writer::AlnWriter;
pub mod idx_check;
#[cfg(not(target_family = "wasm"))]
use crate::ska_ref::idx_check::IdxCheck;

#[cfg(not(target_family = "wasm"))]
use crate::merge_ska_dict::MergeSkaDict;
#[cfg(target_family = "wasm")]
use crate::ska_dict::bit_encoding::UInt;
#[cfg(not(target_family = "wasm"))]
use crate::ska_dict::bit_encoding::{UInt, RC_IUPAC};
use crate::ska_dict::split_kmer::SplitKmer;

#[cfg(target_family = "wasm")]
use crate::logw;
#[cfg(target_family = "wasm")]
use crate::wasm::{fastx_wasm::open_fasta, ska_map::Variant};
#[cfg(target_family = "wasm")]
use seq_io::fasta::Record;
#[cfg(target_family = "wasm")]
use std::io::Read;

/// A split k-mer in the reference sequence encapsulated with positional data.
#[derive(Debug, Clone)]
pub struct RefKmer<IntT> {
    /// Encoded split k-mer
    pub kmer: IntT,
    /// Middle base
    pub base: u8,
    /// Position in the chromosome
    pub pos: usize,
    /// Index of the chromosome
    pub chrom: usize,
    /// Whether on the reverse strand
    pub rc: bool,
}

/// A reference sequence, a list of its split k-mers, and optionally mapping information.
///
/// The default to after building with [`RefSka::new()`] will be a list of split-kmers,
/// as [`Vec<RefKmer>`], along with the reference sequence itself for lookup purposes.
///
/// After running [`RefSka::map()`] against a [`MergeSkaDict`] mapped middle
/// bases and positions will also be populated.
pub struct RefSka<IntT>
where
    IntT: for<'a> UInt<'a>,
{
    /// k-mer size
    k: usize,
    /// Concatenated list of split k-mers
    split_kmer_pos: Vec<RefKmer<IntT>>,
    /// Replace ambiguous bases with N
    ambig_mask: bool,

    /// Input sequence
    /// Chromosome names
    chrom_names: Vec<String>,
    /// Sequence, indexed by chromosome, then position
    seq: Vec<Vec<u8>>,
    /// Where repeats should be masked with 'N'
    repeat_coors: Vec<usize>,

    /// Mapping information
    /// Positions of mapped bases as (chrom, pos)
    #[cfg(not(target_family = "wasm"))]
    mapped_pos: Vec<(usize, usize)>,
    /// Array of mapped bases, rows loci, columns samples
    #[cfg(not(target_family = "wasm"))]
    mapped_variants: Array2<u8>,
    /// Names of the mapped samples
    #[cfg(not(target_family = "wasm"))]
    mapped_names: Vec<String>,
}

/// [`u8`] representation used elsewhere to [`noodles_vcf::record::reference_bases::Base`]
#[cfg(not(target_family = "wasm"))]
#[inline]
fn u8_to_base(ref_base: u8) -> Base {
    match ref_base {
        b'A' => Base::A,
        b'C' => Base::C,
        b'G' => Base::G,
        b'T' => Base::T,
        _ => Base::N,
    }
}

/// The VCF KEYS field used is currently just genotype (GT)
/// These can be used as [`Keys`] in the genotype builder
#[cfg(not(target_family = "wasm"))]
#[inline]
fn gt_keys() -> Keys {
    Keys::try_from(vec![key::GENOTYPE]).unwrap()
}

impl<IntT> RefSka<IntT>
where
    IntT: for<'a> UInt<'a>,
{
    #[cfg(not(target_family = "wasm"))]
    /// Whether [`map`] has been run
    fn is_mapped(&self) -> bool {
        self.mapped_variants.nrows() > 0
    }

    /// Create a list of split k-mers from an input FASTA file.
    ///
    /// Input may have multiple sequences (which are treated as chromosomes and
    /// their indexes maintained).
    ///
    /// # Panics
    /// If an invalid k-mer length (<5, >63 or even) is used.
    ///
    /// Or if input file is invalid:
    /// - File doesn't exist or can't be opened.
    /// - File cannot be parsed as FASTA (FASTQ is not supported).
    /// - If there are no valid split k-mers.
    #[cfg(not(target_family = "wasm"))]
    pub fn new(k: usize, filename: &str, rc: bool, ambig_mask: bool, repeat_mask: bool) -> Self {
        if !(5..=63).contains(&k) || k.is_multiple_of(2) {
            panic!("Invalid k-mer length");
        }

        let mut split_kmer_pos = Vec::new();
        let mut seq = Vec::new();
        let mut chrom_names = Vec::new();
        let mut singles = HashSet::new();
        let mut repeats = HashSet::new();

        let mut reader =
            parse_fastx_file(filename).unwrap_or_else(|_| panic!("Invalid path/file: {filename}",));
        let mut chrom = 0;
        while let Some(record) = reader.next() {
            let seqrec = record.expect("Invalid FASTA record");
            if seqrec.format() == Format::Fastq {
                panic!("Cannot create reference from FASTQ files");
            }
            let chrom_name = str::from_utf8(seqrec.id())
                .unwrap()
                .split_whitespace()
                .next()
                .unwrap();
            chrom_names.push(chrom_name.to_string());
            split_kmer_pos.reserve(seqrec.num_bases());

            let kmer_opt = SplitKmer::new(
                seqrec.seq(),
                seqrec.num_bases(),
                None,
                k,
                rc,
                0,
                QualFilter::NoFilter,
                false,
            );
            if let Some(mut kmer_it) = kmer_opt {
                let (kmer, base, rc) = kmer_it.get_curr_kmer();
                let mut pos = kmer_it.get_middle_pos();
                split_kmer_pos.push(RefKmer {
                    kmer,
                    base,
                    pos,
                    chrom,
                    rc,
                });
                if repeat_mask {
                    Self::track_repeats(kmer, &mut singles, &mut repeats);
                }
                while let Some((kmer, base, rc)) = kmer_it.get_next_kmer() {
                    pos = kmer_it.get_middle_pos();
                    split_kmer_pos.push(RefKmer {
                        kmer,
                        base,
                        pos,
                        chrom,
                        rc,
                    });
                    if repeat_mask {
                        Self::track_repeats(kmer, &mut singles, &mut repeats);
                    }
                }
            }
            chrom += 1;
            seq.push(seqrec.seq().to_vec());
        }
        if split_kmer_pos.is_empty() {
            panic!("{filename} has no valid sequence");
        }

        // Find the repeat ranges, and intersect them
        let mut repeat_coors = Vec::new();
        if repeat_mask {
            let half_split_len = (k - 1) / 2;
            let mut last_chrom = 0;
            let mut last_end = 0;
            let mut chrom_offset = 0;
            for sk in &split_kmer_pos {
                if sk.chrom > last_chrom {
                    chrom_offset += seq[last_chrom].len();
                    last_chrom = sk.chrom;
                }
                if repeats.contains(&sk.kmer) {
                    let start = sk.pos - half_split_len + chrom_offset;
                    let end = sk.pos + half_split_len + chrom_offset;
                    let range = if start > last_end || start == 0 {
                        std::ops::Range {
                            start,
                            end: end + 1,
                        }
                    } else {
                        std::ops::Range {
                            start: last_end + 1,
                            end: end + 1,
                        }
                    };
                    for pos in range {
                        repeat_coors.push(pos);
                    }
                    last_chrom = sk.chrom;
                    last_end = end;
                }
            }
            log::info!(
                "Masking {} unique split k-mer repeats spanning {} bases",
                repeats.len(),
                repeat_coors.len()
            );
        }

        Self {
            k,
            seq,
            ambig_mask,
            chrom_names,
            split_kmer_pos,
            repeat_coors,
            mapped_pos: Vec::new(),
            mapped_variants: Array2::zeros((0, 0)),
            mapped_names: Vec::new(),
        }
    }

    /// Create a list of split k-mers from an input FASTA file.
    ///
    /// Input may have multiple sequences (which are treated as chromosomes and
    /// their indexes maintained).
    ///
    /// # Panics
    /// If an invalid k-mer length (<5, >63 or even) is used.
    ///
    /// Or if input file is invalid:
    /// - File doesn't exist or can't be opened.
    /// - File cannot be parsed as FASTA (FASTQ is not supported).
    /// - If there are no valid split k-mers.
    #[cfg(target_family = "wasm")]
    pub fn new<F: Read>(
        k: usize,
        file: &mut F,
        rc: bool,
        ambig_mask: bool,
        repeat_mask: bool,
    ) -> Self {
        let mut reader = open_fasta(file);

        if !(5..=63).contains(&k) || k.is_multiple_of(2) {
            panic!("Invalid k-mer length");
        }

        let mut split_kmer_pos = Vec::new();
        let mut seq = Vec::new();
        let mut chrom_names = Vec::new();
        let mut singles = HashSet::new();
        let mut repeats = HashSet::new();

        let mut chrom = 0;
        while let Some(record) = reader.next() {
            let seqrec = record.expect("Invalid FASTA record");
            chrom_names.push(seqrec.id().map_or_else(|_| format!("contig_{}", chrom + 1), |id| id.to_string()));
            split_kmer_pos.reserve(
                seqrec
                    .full_seq()
                    .to_vec()
                    .iter()
                    .filter(|&x| *x != 10)
                    .cloned()
                    .collect::<Vec<_>>()
                    .len(),
            );

            let kmer_opt = SplitKmer::new(
                // Remove \n characters from the sequence
                seqrec
                    .seq()
                    .to_vec()
                    .iter()
                    .filter(|&x| *x != 10)
                    .cloned()
                    .collect(),
                seqrec
                    .seq()
                    .to_vec()
                    .iter()
                    .filter(|&x| *x != 10)
                    .cloned()
                    .collect::<Vec<_>>()
                    .len(),
                None,
                k,
                rc,
                0,
                QualFilter::NoFilter,
                false,
            );
            if let Some(mut kmer_it) = kmer_opt {
                let (kmer, base, rc) = kmer_it.get_curr_kmer();
                let mut pos = kmer_it.get_middle_pos();
                split_kmer_pos.push(RefKmer {
                    kmer,
                    base,
                    pos,
                    chrom,
                    rc,
                });
                if repeat_mask {
                    Self::track_repeats(kmer, &mut singles, &mut repeats);
                }
                while let Some((kmer, base, rc)) = kmer_it.get_next_kmer() {
                    pos = kmer_it.get_middle_pos();
                    split_kmer_pos.push(RefKmer {
                        kmer,
                        base,
                        pos,
                        chrom,
                        rc,
                    });
                    if repeat_mask {
                        Self::track_repeats(kmer, &mut singles, &mut repeats);
                    }
                }
            }
            chrom += 1;
            // Remove \n characters from the sequence
            seq.push(
                seqrec
                    .seq()
                    .to_vec()
                    .iter()
                    .filter(|&x| *x != 10)
                    .cloned()
                    .collect::<Vec<_>>(),
            );
        }
        if split_kmer_pos.is_empty() {
            panic!("No valid sequence");
        }

        // Find the repeat ranges, and intersect them
        let mut repeat_coors = Vec::new();
        if repeat_mask {
            let half_split_len = (k - 1) / 2;
            let mut last_chrom = 0;
            let mut last_end = 0;
            let mut chrom_offset = 0;
            for sk in &split_kmer_pos {
                if sk.chrom > last_chrom {
                    chrom_offset += seq[last_chrom].len();
                    last_chrom = sk.chrom;
                }
                if repeats.contains(&sk.kmer) {
                    let start = sk.pos - half_split_len + chrom_offset;
                    let end = sk.pos + half_split_len + chrom_offset;
                    let range = if start > last_end || start == 0 {
                        std::ops::Range {
                            start,
                            end: end + 1,
                        }
                    } else {
                        std::ops::Range {
                            start: last_end + 1,
                            end: end + 1,
                        }
                    };
                    for pos in range {
                        repeat_coors.push(pos);
                    }
                    last_chrom = sk.chrom;
                    last_end = end;
                }
            }
            if cfg!(debug_assertions) {
                logw(
                    &format!(
                        "Masking {} unique split k-mer repeats spanning {} bases",
                        repeats.len(),
                        repeat_coors.len()
                    ),
                    None,
                );
            }
        }

        Self {
            k,
            seq,
            ambig_mask,
            chrom_names,
            split_kmer_pos,
            repeat_coors,
            // mapped_pos: Vec::new(),
            // mapped_variants: Array2::zeros((0, 0)),
            // mapped_names: Vec::new(),
        }
    }

    // Keeps track of split k-mers in the ref, any found before are moved
    // to the repeats set
    fn track_repeats(kmer: IntT, singles: &mut HashSet<IntT>, repeats: &mut HashSet<IntT>) {
        if let Vacant(rep_kmer) = repeats.entry(kmer) {
            match singles.entry(kmer) {
                Vacant(single_entry) => single_entry.insert(),
                Occupied(_) => {
                    rep_kmer.insert();
                }
            }
        }
    }

    /// Map split k-mers in a [`MergeSkaDict`] against the refrence sequence.
    ///
    /// Iterates through the referecne list and searches for split k-mers in
    /// the merged dictionary. These are then added as rows of an array,
    /// with position metadata stored in companion vectors.
    ///
    /// # Panics
    ///
    /// If k-mer sizes are incompatible
    #[cfg(not(target_family = "wasm"))]
    pub fn map(&mut self, ska_dict: &MergeSkaDict<IntT>) {
        if self.k != ska_dict.kmer_len() {
            panic!(
                "K-mer sizes do not match ref:{} skf:{}",
                self.k,
                ska_dict.ksize()
            );
        }
        self.mapped_names = ska_dict.names().clone();
        self.mapped_variants = Array2::zeros((0, ska_dict.nsamples()));
        for ref_k in &self.split_kmer_pos {
            if ska_dict.kmer_dict().contains_key(&ref_k.kmer) {
                let seq_char: Vec<u8> = ska_dict.kmer_dict()[&ref_k.kmer]
                    .iter()
                    .map(|x| match ref_k.rc {
                        true => RC_IUPAC[*x as usize],
                        false => *x,
                    })
                    .collect();
                self.mapped_variants
                    .push_row(ArrayView::from(&seq_char))
                    .unwrap();
                self.mapped_pos.push((ref_k.chrom, ref_k.pos));
            }
        }
    }

    /// Number of valid split k-mers in the reference.
    pub fn ksize(&self) -> usize {
        self.split_kmer_pos.len()
    }

    /// An [`Iterator`] over the reference's split-kmers.
    pub fn kmer_iter(&self) -> impl Iterator<Item = IntT> + '_ {
        self.split_kmer_pos.iter().map(|k| k.kmer)
    }

    /// A direct [`Iterator`] over the reference's split-kmers.
    #[cfg(target_family = "wasm")]
    pub fn kmers(&self) -> impl Iterator<Item = &RefKmer<IntT>> + '_ {
        self.split_kmer_pos.iter()
    }

    #[cfg(not(target_family = "wasm"))]
    /// Calls the necessary parts of AlnWriter (in parallel) to produce all the
    /// pseudoalignments. The calling function either modifies these (VCF) or
    /// simply writes them out (ALN)
    pub fn pseudoalignment(&self, threads: usize) -> Vec<AlnWriter<'_>> {
        if !self.is_mapped() {
            panic!("No split k-mers mapped to reference");
        }
        if self.ambig_mask {
            log::info!("Masking any ambiguous bases (non-A/C/G/T/U/N/-) with 'N'");
        }

        let mut seq_writers =
            vec![
                AlnWriter::new(&self.seq, self.k, &self.repeat_coors, self.ambig_mask);
                self.mapped_names.len()
            ];
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .unwrap();
        seq_writers
            .par_iter_mut()
            .enumerate()
            .for_each(|(idx, seq)| {
                let sample_vars = self.mapped_variants.slice(s![.., idx]);
                for ((mapped_chrom, mapped_pos), base) in
                    self.mapped_pos.iter().zip(sample_vars.iter())
                {
                    if *base != b'-' {
                        seq.write_split_kmer(*mapped_pos, *mapped_chrom, *base);
                    }
                }
                seq.finalise();
            });
        seq_writers
    }

    #[cfg(target_family = "wasm")]
    /// Calls the necessary parts of AlnWriter (in parallel) to produce all the
    /// pseudoalignments. The calling function simply writes them out (ALN)
    pub fn pseudoalignment(&self, mapped_bases: &[Variant]) -> Vec<String> {
        let mapped_variants: Vec<u8> = mapped_bases.iter().map(|v| v.base).collect();
        let mapped_pos: Vec<(usize, usize)> =
            mapped_bases.iter().map(|v| (v.chrom, v.pos)).collect();

        let mut seq_writers =
            vec![AlnWriter::new(&self.seq, self.k, &self.repeat_coors, self.ambig_mask); 1];
        seq_writers
            .par_iter_mut()
            .enumerate()
            .for_each(|(_idx, seq)| {
                let sample_vars = mapped_variants.clone();
                for ((mapped_chrom, mapped_pos), base) in mapped_pos.iter().zip(sample_vars.iter())
                {
                    if *base != b'-' {
                        seq.write_split_kmer(*mapped_pos, *mapped_chrom, *base);
                    }
                }
                seq.finalise()
            });

        let sequences: Vec<String> = seq_writers
            .iter_mut()
            .map(|seq| String::from_utf8_lossy(seq.get_seq()).to_string())
            .collect();

        sequences
    }

    /// Write mapped variants as an aln/multi-FASTA file, using [`needletail`].
    ///
    /// Rows are samples, columns are variants, so this is broadly equivalent
    /// to writing out the transpose of the mapped variant array.
    ///
    /// Extra work is used to check to write out flanking bases of the split
    /// k-mer, and padding with any missing positions.
    ///
    /// This is done in memory, then samples are written to the output buffer
    /// one at a time.
    ///
    /// # Panics
    ///
    /// If [`RefSka::map()`] has not been run yet, or no split-kmers mapped.
    #[cfg(not(target_family = "wasm"))]
    pub fn write_aln<W: Write>(
        &self,
        f: &mut W,
        threads: usize,
    ) -> Result<(), needletail::errors::ParseError> {
        if self.chrom_names.len() > 1 {
            log::warn!(
                "Reference contained multiple contigs, in the output they will be concatenated"
            );
        }

        let alignments = self.pseudoalignment(threads);
        log::info!("Writing alignment to file");
        for (sample_name, mut aligned_seq) in self.mapped_names.iter().zip(alignments) {
            write_fasta(
                sample_name.as_bytes(),
                aligned_seq.get_seq(),
                f,
                needletail::parser::LineEnding::Unix,
            )?;
        }
        Ok(())
    }

    /// Write mapped variants as a VCF, using [`noodles_vcf`].
    ///
    /// This uses [`RefSka::write_aln()`] and converts the output to a VCF:
    /// - A basic VCF header is written out first.
    /// - An alignment is created (in memory), which is then converted to a VCF
    /// - Any ambiguous base is currently output as N, for simplicity of GT field.
    ///
    /// # Panics
    ///
    /// If [`RefSka::map()`] has not been run yet, or no split-kmers mapped to
    /// the reference.
    #[cfg(not(target_family = "wasm"))]
    pub fn write_vcf<W: Write>(&self, f: &mut W, threads: usize) -> Result<(), std::io::Error> {
        if !self.is_mapped() {
            panic!("No split k-mers mapped to reference");
        }

        // Get the pseudoalignment and store in array form
        let alignments = self.pseudoalignment(threads);

        log::info!("Converting to VCF");
        let mut variants = Array2::zeros((0, alignments[0].total_size()));
        for mut seq in alignments {
            variants.push_row(ArrayView::from(seq.get_seq())).unwrap();
        }
        // Transpose the array
        let var_t = variants.t();
        let mut var_t_owned = Array2::zeros(var_t.raw_dim());
        var_t_owned.assign(&var_t);

        // Write the VCF header
        let mut writer = vcf::Writer::new(f);
        let mut header_builder = vcf::Header::builder();
        for contig in &self.chrom_names {
            header_builder = header_builder.add_contig(
                contig.parse().expect("Could not add contig to header"),
                Map::<Contig>::new(),
            );
        }
        for name in &self.mapped_names {
            header_builder = header_builder.add_sample_name(name);
        }
        let header = header_builder.build();
        writer.write_header(&header)?;

        // Write each record (column)
        let keys = gt_keys();
        let idx_map = IdxCheck::new(&self.seq);
        for (sample_variants, (map_chrom, map_pos)) in var_t_owned.outer_iter().zip(idx_map.iter())
        {
            let ref_base = self.seq[map_chrom][map_pos];
            let ref_allele = u8_to_base(ref_base);

            let mut genotype_vec = Vec::with_capacity(var_t_owned.ncols());
            let mut alt_bases: Vec<Base> = Vec::new();
            let mut variant = false;
            for mapped_base in sample_variants {
                let gt = if *mapped_base == ref_base {
                    String::from("0")
                } else if *mapped_base == b'-' {
                    variant = true;
                    String::from(".")
                } else {
                    variant = true;
                    let alt_base = u8_to_base(*mapped_base);
                    if !alt_bases.contains(&alt_base) {
                        alt_bases.push(alt_base);
                    }
                    (alt_bases.iter().position(|&r| r == alt_base).unwrap() + 1).to_string()
                };
                genotype_vec.push(vec![Some(Value::String(gt))]);
            }
            if variant {
                let genotypes = Genotypes::new(keys.clone(), genotype_vec);
                let alt_alleles: Vec<Allele> =
                    alt_bases.iter().map(|a| Allele::Bases(vec![*a])).collect();
                let record = vcf::Record::builder()
                    .set_chromosome(
                        self.chrom_names[map_chrom]
                            .parse()
                            .expect("Invalid chromosome name"),
                    )
                    .set_position(Position::from(map_pos + 1))
                    .add_reference_base(ref_allele)
                    .set_alternate_bases(AlternateBases::from(alt_alleles))
                    .set_genotypes(genotypes)
                    .build()
                    .expect("Could not construct record");
                writer.write_record(&header, &record)?;
            }
        }
        Ok(())
    }

    #[cfg(target_family = "wasm")]
    fn sanitise_vcf_sample_name(name: &str) -> String {
        name.chars()
            .map(|c| match c {
                '\t' | '\r' | '\n' => '_',
                _ => c,
            })
            .collect()
    }

    #[cfg(target_family = "wasm")]
    fn normalise_vcf_base(base: u8) -> u8 {
        match base.to_ascii_uppercase() {
            b'A' | b'C' | b'G' | b'T' => base.to_ascii_uppercase(),
            _ => b'N',
        }
    }

    #[cfg(target_family = "wasm")]
    fn format_wasm_gt(base: u8, ref_base: u8, alt_bases: &mut Vec<u8>) -> String {
        if base == b'-' {
            return String::from(".");
        }

        let normalised_base = Self::normalise_vcf_base(base);
        if normalised_base == ref_base {
            return String::from("0");
        }

        if let Some(idx) = alt_bases.iter().position(|alt| *alt == normalised_base) {
            (idx + 1).to_string()
        } else {
            alt_bases.push(normalised_base);
            alt_bases.len().to_string()
        }
    }

    #[cfg(target_family = "wasm")]
    /// Writes mapped variants for browser use as combined VCF text.
    pub fn write_vcf_text_from_sparse_maps(
        &self,
        sample_names: &[String],
        sample_maps: &[Vec<Variant>],
    ) -> String {
        let mut mapped_by_sample: Vec<HashMap<(usize, usize), u8>> = Vec::with_capacity(sample_maps.len());
        for sample_map in sample_maps {
            let mut mapped = HashMap::with_capacity(sample_map.len());
            for variant in sample_map {
                mapped.insert((variant.chrom, variant.pos), variant.base);
            }
            mapped_by_sample.push(mapped);
        }

        let mut vcf = String::new();
        vcf.push_str("##fileformat=VCFv4.3\n");
        vcf.push_str("##source=ska.rust\n");
        for chrom_name in &self.chrom_names {
            vcf.push_str("##contig=<ID=");
            vcf.push_str(chrom_name);
            vcf.push_str(">\n");
        }
        vcf.push_str("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
        vcf.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
        for sample_name in sample_names {
            vcf.push('\t');
            vcf.push_str(&Self::sanitise_vcf_sample_name(sample_name));
        }
        vcf.push('\n');

        for ref_kmer in self.kmers() {
            let ref_base = Self::normalise_vcf_base(self.seq[ref_kmer.chrom][ref_kmer.pos]);
            let mut alt_bases: Vec<u8> = Vec::new();
            let mut genotypes: Vec<String> = Vec::with_capacity(mapped_by_sample.len());

            for mapped in &mapped_by_sample {
                let base = mapped
                    .get(&(ref_kmer.chrom, ref_kmer.pos))
                    .copied()
                    .unwrap_or(b'-');
                genotypes.push(Self::format_wasm_gt(base, ref_base, &mut alt_bases));
            }

            if alt_bases.is_empty() {
                continue;
            }

            let chrom_name = self
                .chrom_names
                .get(ref_kmer.chrom)
                .map(String::as_str)
                .unwrap_or("contig");
            vcf.push_str(chrom_name);
            vcf.push('\t');
            vcf.push_str(&(ref_kmer.pos + 1).to_string());
            vcf.push_str("\t.\t");
            vcf.push(ref_base as char);
            vcf.push('\t');
            for (idx, alt_base) in alt_bases.iter().enumerate() {
                if idx > 0 {
                    vcf.push(',');
                }
                vcf.push(*alt_base as char);
            }
            vcf.push_str("\t.\t.\t.\tGT");
            for genotype in genotypes {
                vcf.push('\t');
                vcf.push_str(&genotype);
            }
            vcf.push('\n');
        }

        vcf
    }

    /// Returns the k value associated with this SkaRef struct
    #[cfg(target_family = "wasm")]
    pub fn get_k(&self) -> usize {
        self.k
    }

    /// Returns a reference to the sequence
    #[cfg(target_family = "wasm")]
    pub fn get_seq(&self) -> &Vec<Vec<u8>> {
        &self.seq
    }
}
