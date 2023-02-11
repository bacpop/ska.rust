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
//! let ska_dict = load_array(&["tests/test_files_in/merge.skf".to_string()], threads).to_dict();
//!
//! // Index a reference sequence
//! let mut ref_kmers = RefSka::new(ska_dict.kmer_len(), &"tests/test_files_in/test_ref.fa", ska_dict.rc());
//!
//! // Run mapping, output an alignment to stdout
//! ref_kmers.map(&ska_dict);
//! let mut out_stream = set_ostream(&None);
//! ref_kmers.write_aln(&mut out_stream, threads);
//! ```

use std::io::Write;
use std::str;

use rayon::prelude::*;

use noodles_vcf::{
    self as vcf,
    header::format::Key,
    header::record::value::{map::Contig, Map},
    record::{
        alternate_bases::Allele,
        genotypes::{
            genotype::{field::Value, Field},
            Genotype, Keys,
        },
        reference_bases::Base,
        AlternateBases, Genotypes, Position,
    },
};

extern crate needletail;
use ndarray::{s, Array2, ArrayView};
use needletail::{
    parse_fastx_file,
    parser::{write_fasta, Format},
};

pub mod aln_writer;
use crate::ska_ref::aln_writer::AlnWriter;
pub mod idx_check;
use crate::ska_ref::idx_check::IdxCheck;

use crate::merge_ska_dict::MergeSkaDict;
use crate::ska_dict::bit_encoding::{RC_IUPAC, RevComp};
use crate::ska_dict::split_kmer::SplitKmer;

/// A split k-mer in the reference sequence encapsulated with positional data.
#[derive(Debug, Clone)]
pub struct RefKmer<IntT>
where
    IntT: RevComp
{
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
    IntT: RevComp
{
    /// k-mer size
    k: usize,
    /// Concatenated list of split k-mers
    split_kmer_pos: Vec<RefKmer<IntT>>,

    /// Input sequence

    /// Chromosome names
    chrom_names: Vec<String>,
    /// Sequence, indexed by chromosome, then position
    seq: Vec<Vec<u8>>,

    /// Mapping information

    /// Positions of mapped bases as (chrom, pos)
    mapped_pos: Vec<(usize, usize)>,
    /// Array of mapped bases, rows loci, columns samples
    mapped_variants: Array2<u8>,
    /// Names of the mapped samples
    mapped_names: Vec<String>,
}

/// [`u8`] representation used elsewhere to [`noodles_vcf::record::reference_bases::Base`]
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
#[inline]
fn gt_keys() -> Keys {
    "GT".parse().expect("Genotype format error")
}

impl<IntT> RefSka<IntT>
where
    IntT: RevComp
{
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
    /// If an invalid k-mer length (<5, >31 or even) is used.
    ///
    /// Or if input file is invalid:
    /// - File doesn't exist or can't be opened.
    /// - File cannot be parsed as FASTA (FASTQ is not supported).
    /// - If there are no valid split k-mers.
    pub fn new(k: usize, filename: &str, rc: bool) -> Self {
        if !(5..=31).contains(&k) || k % 2 == 0 {
            panic!("Invalid k-mer length");
        }

        let mut split_kmer_pos = Vec::new();
        let mut seq = Vec::new();
        let mut chrom_names = Vec::new();

        let mut reader =
            parse_fastx_file(filename).unwrap_or_else(|_| panic!("Invalid path/file: {filename}",));
        let mut chrom = 0;
        while let Some(record) = reader.next() {
            let seqrec = record.expect("Invalid FASTA record");
            if seqrec.format() == Format::Fastq {
                panic!("Cannot create reference from FASTQ files");
            }
            chrom_names.push(str::from_utf8(seqrec.id()).unwrap().to_owned());
            split_kmer_pos.reserve(seqrec.num_bases());

            let kmer_opt = SplitKmer::new(seqrec.seq(), seqrec.num_bases(), None, k, rc, 0);
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
                while let Some((kmer, base, rc)) = kmer_it.get_next_kmer() {
                    pos = kmer_it.get_middle_pos();
                    split_kmer_pos.push(RefKmer {
                        kmer,
                        base,
                        pos,
                        chrom,
                        rc,
                    });
                }
            }
            chrom += 1;
            seq.push(seqrec.seq().to_vec());
        }
        if split_kmer_pos.is_empty() {
            panic!("{filename} has no valid sequence");
        }

        Self {
            k,
            seq,
            chrom_names,
            split_kmer_pos,
            mapped_pos: Vec::new(),
            mapped_variants: Array2::zeros((0, 0)),
            mapped_names: Vec::new(),
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

    // Calls the necessary parts of AlnWriter (in parallel) to produce all the
    // pseudoalignments. The calling function either modifies these (VCF) or
    // simply writes them out (ALN)
    fn pseudoalignment(&self, threads: usize) -> Vec<AlnWriter> {
        if !self.is_mapped() {
            panic!("No split k-mers mapped to reference");
        }

        let mut seq_writers = vec![AlnWriter::new(&self.seq, self.k); self.mapped_names.len()];
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
    pub fn write_aln<W: Write>(
        &self,
        f: &mut W,
        threads: usize,
    ) -> Result<(), needletail::errors::ParseError> {
        if self.chrom_names.len() > 1 {
            eprintln!("WARNING: Reference contained multiple contigs, in the output they will be concatenated");
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
            header_builder = header_builder.add_contig(Map::<Contig>::new(
                contig.parse().expect("Could not add contig to header"),
            ));
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

            let mut genotype_vec = Vec::new();
            genotype_vec.reserve(var_t_owned.ncols());
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
                let field = Field::new(Key::Genotype, Some(Value::String(gt)));
                genotype_vec
                    .push(Genotype::try_from(vec![field]).expect("Could not construct genotypes"));
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
                writer.write_record(&record)?;
            }
        }
        Ok(())
    }
}
