// Class for split-kmers from one input file
// Used for ska map

use std::str;

use std::io::Write;
use noodles_vcf::{
    self as vcf,
    header::record::value::{map::Contig, Map},
    header::format::Key,
    record::{
        reference_bases::Base,
        alternate_bases::Allele, AlternateBases,
        genotypes::{genotype::{field::Value, Field}, Keys, Genotype},
        Genotypes, Position,
    },
};

extern crate needletail;
use needletail::{parse_fastx_file, parser::write_fasta};
use ndarray::{ArrayView, Array2, s};

use crate::merge_ska_dict::MergeSkaDict;
use crate::ska_dict::split_kmer::SplitKmer;
use crate::ska_dict::bit_encoding::RC_IUPAC;

pub struct RefKmer {
    pub kmer: u64,
    pub base: u8,
    pub pos: usize,
    pub chrom: usize,
    pub rc: bool,
}

pub struct RefSka {
    // Index
    k: usize,
    split_kmer_pos: Vec<RefKmer>,

    // Sequence
    chrom_names: Vec<String>,
    seq: Vec<Vec<u8>>,
    total_size: usize,

    // Mapping
    mapped_pos: Vec<(usize, usize)>, // (chrom, pos)
    mapped_variants: Array2<u8>,
    mapped_names: Vec<String>
}

impl RefSka {
    fn is_mapped(&self) -> bool {
        self.mapped_variants.nrows() > 0
    }

    pub fn ksize(&self) -> usize {
        self.split_kmer_pos.len()
    }

    pub fn new(k: usize, filename: &str, rc: bool) -> Self {
        if k < 5 || k > 31 || k % 2 == 0 {
            panic!("Invalid k-mer length");
        }

        let mut split_kmer_pos = Vec::new();
        let mut seq = Vec::new();
        let mut chrom_names = Vec::new();
        let mut total_size = 0;

        let mut reader =
            parse_fastx_file(&filename).expect(&format!("Invalid path/file: {}", filename));
        let mut chrom = 0;
        while let Some(record) = reader.next() {
            let seqrec = record.expect("Invalid FASTA record");
            chrom_names.push(str::from_utf8(seqrec.id()).unwrap().to_owned());
            split_kmer_pos.reserve(seqrec.num_bases());
            total_size += seqrec.num_bases();

            let kmer_opt = SplitKmer::new(seqrec.seq(), seqrec.num_bases(), k, rc);
            if kmer_opt.is_some() {
                let mut kmer_it = kmer_opt.unwrap();
                let (kmer, base, rc) = kmer_it.get_curr_kmer();
                let mut pos = kmer_it.get_pos();
                split_kmer_pos.push(RefKmer{kmer, base, pos, chrom, rc});
                while let Some((kmer, base, rc)) = kmer_it.get_next_kmer() {
                    pos = kmer_it.get_pos();
                    split_kmer_pos.push(RefKmer{kmer, base, pos, chrom, rc});
                }
            }
            chrom += 1;
            seq.push(seqrec.seq().to_vec());
        }
        if split_kmer_pos.len() == 0 {
            panic!("{} has no valid sequence", filename);
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
            mapped_names
        }
    }

    pub fn map(&mut self, ska_dict: &MergeSkaDict) {
        if self.k != ska_dict.ksize() {
            panic!("K-mer sizes do not match ref:{} skf:{}", self.k, ska_dict.ksize());
        }
        self.mapped_names = ska_dict.names().clone();
        self.mapped_variants = Array2::zeros((0, ska_dict.nsamples()));
        for ref_k in &self.split_kmer_pos {
            if ska_dict.kmer_dict().contains_key(&ref_k.kmer) {
                let seq_char: Vec::<u8> = ska_dict.kmer_dict()[&ref_k.kmer]
                    .iter()
                    .map(|x| {
                        match ref_k.rc {
                            true => RC_IUPAC[*x as usize],
                            false => *x
                        }
                    })
                    .collect();
                self.mapped_variants.push_row(ArrayView::from(&seq_char)).unwrap();
                self.mapped_pos.push((ref_k.chrom, ref_k.pos));
            }
        }
    }

    pub fn write_vcf<W: Write>(&self, f: &mut W) -> Result<(), std::io::Error> {
        if !self.is_mapped() {
            panic!("Tried to write VCF before variants mapped");
        }

        // Write the VCF header
        let mut writer = vcf::Writer::new(f);
        let mut header_builder = vcf::Header::builder();
        for contig in &self.chrom_names {
            header_builder = header_builder.add_contig(Map::<Contig>::new(contig.parse().expect("Could not add contig to header")));
        }
        for name in &self.mapped_names {
            header_builder = header_builder.add_sample_name(name);
        }
        let header = header_builder.build();
        writer.write_header(&header)?;

        // Write each record (column)
        let keys: Keys = "GT".parse().expect("Genotype format error");
        for ((map_chrom, map_pos), bases) in self.mapped_pos.iter().zip(self.mapped_variants.outer_iter()) {
            let ref_base = self.seq[*map_chrom][*map_pos];
            let ref_allele = match ref_base {
                b'A' => Base::A,
                b'C' => Base::C,
                b'G' => Base::G,
                b'T' => Base::T,
                _ => Base::N
            };

            let mut genotype_vec = Vec::new();
            genotype_vec.reserve(bases.len());
            let mut alt_bases: Vec<Base> = Vec::new();
            for mapped_base in bases {
                let mut gt: usize = 0;
                if *mapped_base != ref_base {
                    gt = match *mapped_base {
                        b'A' => {
                            if !alt_bases.contains(&Base::A) {
                                alt_bases.push(Base::A);
                            }
                            alt_bases.iter().position(|&r| r == Base::A).unwrap()
                        },
                        b'C' => {
                            if !alt_bases.contains(&Base::C) {
                                alt_bases.push(Base::C);
                            }
                            alt_bases.iter().position(|&r| r == Base::C).unwrap()
                        },
                        b'G' => {
                            if !alt_bases.contains(&Base::G) {
                                alt_bases.push(Base::G);
                            }
                            alt_bases.iter().position(|&r| r == Base::G).unwrap()
                        },
                        b'T' => {
                            if !alt_bases.contains(&Base::T) {
                                alt_bases.push(Base::T);
                            }
                            alt_bases.iter().position(|&r| r == Base::T).unwrap()
                        },
                        _ => {
                            if !alt_bases.contains(&Base::N) {
                                alt_bases.push(Base::N);
                            }
                            alt_bases.iter().position(|&r| r == Base::N).unwrap()
                        },
                    }
                }
                let field = Field::new(Key::Genotype, Some(Value::String(gt.to_string())));
                genotype_vec.push(field);
            }
            let genotypes = Genotypes::new(
                keys.clone(),
                vec![Genotype::try_from(genotype_vec).expect("Could not construct genotypes")],
            );
            let record = vcf::Record::builder()
                .set_chromosome(self.chrom_names[*map_chrom].parse().expect("Invalid chromosome name"))
                .set_position(Position::from(*map_pos))
                .add_reference_base(ref_allele)
                .set_alternate_bases(AlternateBases::from(vec![Allele::Bases(alt_bases)]))
                .set_genotypes(genotypes)
                .build().expect("Could not construct record");
            writer.write_record(&record)?;
        }
        Ok(())
    }

    pub fn write_aln<W: Write>(&self, f: &mut W) -> Result<(), needletail::errors::ParseError> {
        if !self.is_mapped() {
            panic!("Tried to write VCF before variants mapped");
        }
        if self.chrom_names.len() > 1 {
            eprintln!("WARNING: Reference contained multiple contigs, in the output they will be concatenated");
        }
        for (sample_idx, sample_name) in self.mapped_names.iter().enumerate() {
            let sample_vars = self.mapped_variants.slice(s![.., sample_idx]);
            let mut seq: Vec<u8> = Vec::new();
            seq.reserve(self.total_size);

            let (mut next_pos, mut curr_chrom) = (0, 0);
            for ((map_chrom, map_pos), base) in self.mapped_pos.iter().zip(sample_vars.iter()) {
                // Move forward to next chromosome/contig
                if *map_chrom > curr_chrom {
                    seq.extend_from_slice(&self.seq[curr_chrom][next_pos..]);
                    curr_chrom += 1;
                    next_pos = 0;
                }
                if *map_pos > next_pos {
                    // Copy in ref seq if no k-mers mapped over a region
                    seq.extend_from_slice(&self.seq[curr_chrom][next_pos..*map_pos]);
                }
                next_pos = *map_pos + 1;
                seq.push(*base);
            }
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
