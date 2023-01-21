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
//! ref_kmers.write_aln(&mut out_stream);
//! ```

use std::str;

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
    Record,
};
use std::io::Write;

extern crate needletail;
use ndarray::{s, Array2, ArrayView};
use needletail::{
    parse_fastx_file,
    parser::{write_fasta, Format},
};

use crate::merge_ska_dict::MergeSkaDict;
use crate::ska_dict::bit_encoding::RC_IUPAC;
use crate::ska_dict::split_kmer::SplitKmer;

/// A split k-mer in the reference sequence encapsulated with positional data.
#[derive(Debug, Clone)]
pub struct RefKmer {
    /// Encoded split k-mer
    pub kmer: u64,
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
pub struct RefSka {
    /// k-mer size
    k: usize,
    /// Concatenated list of split k-mers
    split_kmer_pos: Vec<RefKmer>,

    /// Input sequence

    /// Chromosome names
    chrom_names: Vec<String>,
    /// Sequence, indexed by chromosome, then position
    seq: Vec<Vec<u8>>,
    /// Total length of all input sequence
    total_size: usize,

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

/// Genotypes "." when no reference mapped to (of sample length)
/// These can be used as genotypes in the record builder
fn generate_missing_genotypes(n_samples: usize) -> Genotypes {
    let missing_field = Field::new(Key::Genotype, Some(Value::String(".".to_string())));
    let missing_genotype_vec = vec![
        Genotype::try_from(vec![missing_field])
            .expect("Could not construct genotypes");
        n_samples
    ];
    Genotypes::new(gt_keys(), missing_genotype_vec)
}

impl RefSka {
    /// Whether [`map`] has been run
    fn is_mapped(&self) -> bool {
        self.mapped_variants.nrows() > 0
    }

    /// An [`Iterator`] between start and end positions which generates records with the passed genotypes
    ///
    /// The result can then be used with a writer. So far, only over missing records.
    ///
    /// (I wrote it as an iterator to avoid passing the writer outside of the write function
    /// some pretty new rust stuff here! ...lifetimes, move, iterator trait)
    fn iter_missing_vcf_rows<'a>(
        &'a self,
        chrom: usize,
        start: usize,
        end: usize,
        geno: &'a Genotypes,
    ) -> impl Iterator<Item = Record> + 'a {
        (start..end).into_iter().map(move |missing_pos| {
            let ref_allele = u8_to_base(self.seq[chrom][missing_pos]);
            vcf::Record::builder()
                .set_chromosome(
                    self.chrom_names[chrom]
                        .parse()
                        .expect("Invalid chromosome name"),
                )
                .set_position(Position::from(missing_pos + 1))
                .add_reference_base(ref_allele)
                .set_genotypes(geno.clone())
                .build()
                .expect("Could not construct record")
        })
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
        let mut total_size = 0;

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
            total_size += seqrec.num_bases();

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

        let mapped_variants = Array2::zeros((0, 0));
        let mapped_pos = Vec::new();
        let mapped_names = Vec::new();
        Self {
            k,
            seq,
            total_size,
            chrom_names,
            split_kmer_pos,
            mapped_pos,
            mapped_variants,
            mapped_names,
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
    pub fn map(&mut self, ska_dict: &MergeSkaDict) {
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
    pub fn kmer_iter(&self) -> impl Iterator<Item = u64> + '_ {
        self.split_kmer_pos.iter().map(|k| k.kmer)
    }

    /// Write mapped variants as a VCF, using [`noodles_vcf`].
    ///
    /// Rows are variants, columns are samples, so this is broadly equivalent
    /// to writing out the mapped variant array and its metadata.
    ///
    /// Extra work is used to check for unmapped missing bases between split k-mer
    /// matches.
    ///
    /// A basic VCF header is written out first.
    ///
    /// Variants are then written to the output buffer one at a time.
    ///
    /// # Panics
    ///
    /// If [`RefSka::map()`] has not been run yet, or no split-kmers mapped
    pub fn write_vcf<W: Write>(&self, f: &mut W) -> Result<(), std::io::Error> {
        if !self.is_mapped() {
            panic!("No split k-mers mapped to reference");
        }

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
        let missing_genotypes = generate_missing_genotypes(self.mapped_names.len());

        // Note that a lot of the logic here is similar to write_aln below, which
        // is a bit simpler to follow
        let half_split_len = (self.k - 1) / 2;
        let (mut next_pos, mut curr_chrom) = (0, 0);
        for ((map_chrom, map_pos), bases) in self
            .mapped_pos
            .iter()
            .zip(self.mapped_variants.outer_iter())
        {
            // Fill missing bases to the end of the last chromosome
            if *map_chrom > curr_chrom {
                for record in self.iter_missing_vcf_rows(
                    curr_chrom,
                    next_pos + half_split_len,
                    self.seq[curr_chrom].len(),
                    &missing_genotypes,
                ) {
                    writer.write_record(&record)?;
                }
                curr_chrom += 1;
                next_pos = 0;
            }
            // Fill missing bases between k-mer matches
            if *map_pos > next_pos {
                for record in self.iter_missing_vcf_rows(
                    *map_chrom,
                    next_pos,
                    *map_pos - half_split_len,
                    &missing_genotypes,
                ) {
                    writer.write_record(&record)?;
                }
            }

            let ref_base = self.seq[*map_chrom][*map_pos];
            let ref_allele = u8_to_base(ref_base);

            let mut genotype_vec = Vec::new();
            genotype_vec.reserve(bases.len());
            let mut alt_bases: Vec<Base> = Vec::new();
            for mapped_base in bases {
                let mut gt = String::from("0");
                if *mapped_base != ref_base {
                    if *mapped_base == b'-' {
                        gt = ".".to_string();
                    } else {
                        let alt_base = u8_to_base(*mapped_base);
                        if !alt_bases.contains(&alt_base) {
                            alt_bases.push(alt_base);
                        }
                        gt = (alt_bases.iter().position(|&r| r == alt_base).unwrap() + 1)
                            .to_string();
                    }
                }
                let field = Field::new(Key::Genotype, Some(Value::String(gt)));
                genotype_vec
                    .push(Genotype::try_from(vec![field]).expect("Could not construct genotypes"));
            }
            if !alt_bases.is_empty() {
                let genotypes = Genotypes::new(keys.clone(), genotype_vec);
                let alt_alleles: Vec<Allele> =
                    alt_bases.iter().map(|a| Allele::Bases(vec![*a])).collect();
                let record = vcf::Record::builder()
                    .set_chromosome(
                        self.chrom_names[*map_chrom]
                            .parse()
                            .expect("Invalid chromosome name"),
                    )
                    .set_position(Position::from(*map_pos + 1))
                    .add_reference_base(ref_allele)
                    .set_alternate_bases(AlternateBases::from(alt_alleles))
                    .set_genotypes(genotypes)
                    .build()
                    .expect("Could not construct record");
                writer.write_record(&record)?;
            }
            next_pos = *map_pos + 1;
        }
        // Fill any missing bases at the end of final contig
        let final_chrom_id = self.chrom_names.len() - 1;
        for record in self.iter_missing_vcf_rows(
            final_chrom_id,
            next_pos + half_split_len,
            self.seq[final_chrom_id].len(),
            &missing_genotypes,
        ) {
            writer.write_record(&record)?;
        }
        Ok(())
    }

    /// Write mapped variants as an aln/multi-FASTA file, using [`needletail`].
    ///
    /// Rows are samples, columns are variants, so this is broadly equivalent
    /// to writing out the transpose of the mapped variant array.
    ///
    /// Extra work is used to check to write out flanking bases of the split
    /// k-mer, and padding with any missing positions.
    ///
    /// Samples are written to the output buffer one at a time.
    ///
    /// # Panics
    ///
    /// If [`RefSka::map()`] has not been run yet, or no split-kmers mapped
    pub fn write_aln<W: Write>(&self, f: &mut W) -> Result<(), needletail::errors::ParseError> {
        if !self.is_mapped() {
            panic!("TNo split k-mers mapped to reference");
        }
        if self.chrom_names.len() > 1 {
            eprintln!("WARNING: Reference contained multiple contigs, in the output they will be concatenated");
        }
        let half_split_len = (self.k - 1) / 2;
        for (sample_idx, sample_name) in self.mapped_names.iter().enumerate() {
            let sample_vars = self.mapped_variants.slice(s![.., sample_idx]);
            let mut seq: Vec<u8> = Vec::new();
            seq.reserve(self.total_size);

            // if this proves tricky, it might be better to iterate through non-missing
            // matches and paste each flanking end and middle base into the right place
            // in the vec (and skip to next match where right end is beyond last paste)
            // 	the alternative would be to:
            //      - allocate all missing
            //      - iterate through ref
            //      - write the whole k-mer in when found (ref, base, ref)
            //      - then skip over k in ref
            // (even just modifying map to skip over k would work, but not a VCF)
            let (mut next_pos, mut curr_chrom) = (0, 0);
            for ((map_chrom, map_pos), base) in self.mapped_pos.iter().zip(sample_vars.iter()) {
                // Move forward to next chromosome/contig
                if *map_chrom > curr_chrom {
                    seq.extend_from_slice(
                        &self.seq[curr_chrom][next_pos..(next_pos + half_split_len)],
                    );
                    seq.extend(vec![
                        b'-';
                        self.seq[curr_chrom].len() - (next_pos + half_split_len)
                    ]);
                    curr_chrom += 1;
                    next_pos = 0;
                }
                if *map_pos > next_pos {
                    // Missing bases, if no k-mers mapped over a region
                    seq.extend(vec![b'-'; *map_pos - next_pos - half_split_len]);
                    seq.extend_from_slice(
                        &self.seq[curr_chrom][(*map_pos - half_split_len)..*map_pos],
                    );
                }
                next_pos = *map_pos + 1;
                if *base == b'-' {
                    // This is around a split k-mer match, so we can fill in the
                    // flanking region with the reference
                    seq.push(self.seq[curr_chrom][*map_pos]);
                } else {
                    seq.push(*base);
                }
            }
            // Fill up to end of contig
            seq.extend_from_slice(&self.seq[curr_chrom][next_pos..(next_pos + half_split_len)]);
            seq.extend(vec![
                b'-';
                self.seq[curr_chrom].len() - (next_pos + half_split_len)
            ]);
            write_fasta(
                sample_name.as_bytes(),
                seq.as_slice(),
                f,
                needletail::parser::LineEnding::Unix,
            )?;
        }
        Ok(())
    }
}
