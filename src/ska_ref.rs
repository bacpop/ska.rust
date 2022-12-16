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

fn u8_to_Base(ref_base: u8) -> Base {
    match ref_base {
        b'A' => Base::A,
        b'C' => Base::C,
        b'G' => Base::G,
        b'T' => Base::T,
        _ => Base::N
    }
}

impl RefSka {
    fn is_mapped(&self) -> bool {
        self.mapped_variants.nrows() > 0
    }

    pub fn ksize(&self) -> usize {
        self.split_kmer_pos.len()
    }

    pub fn kmer_iter(&self) -> impl Iterator<Item=u64> + '_ {
        self.split_kmer_pos.iter().map(|k| k.kmer)
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
                let mut pos = kmer_it.get_middle_pos();
                split_kmer_pos.push(RefKmer{kmer, base, pos, chrom, rc});
                while let Some((kmer, base, rc)) = kmer_it.get_next_kmer() {
                    pos = kmer_it.get_middle_pos();
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
        if self.k != ska_dict.kmer_len() {
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
        // TODO: need to add in missing positions
        // TODO: deal with missing positions over chromosomes and at ends
        let keys: Keys = "GT".parse().expect("Genotype format error");
        let missing_field = Field::new(Key::Genotype, Some(Value::String(".".to_string())));
        let missing_genotype_vec = vec![Genotype::try_from(vec![missing_field]).expect("Could not construct genotypes"); self.mapped_names.len()];
        let missing_genotypes = Genotypes::new(
            keys.clone(),
            missing_genotype_vec,
        );

        let (mut next_pos, mut curr_chrom) = (0, 0);
        for ((map_chrom, map_pos), bases) in self.mapped_pos.iter().zip(self.mapped_variants.outer_iter()) {
            if *map_pos > next_pos {
                for missing_pos in next_pos..*map_pos {
                    let ref_allele = u8_to_Base(self.seq[*map_chrom][missing_pos]);
                    let record = vcf::Record::builder()
                        .set_chromosome(self.chrom_names[*map_chrom].parse().expect("Invalid chromosome name"))
                        .set_position(Position::from(missing_pos))
                        .add_reference_base(ref_allele)
                        .set_genotypes(missing_genotypes.clone())
                        .build().expect("Could not construct record");
                    writer.write_record(&record)?;
                }
            }

            let ref_base = self.seq[*map_chrom][*map_pos];
            let ref_allele = u8_to_Base(ref_base);

            let mut genotype_vec = Vec::new();
            genotype_vec.reserve(bases.len());
            let mut alt_bases: Vec<Base> = Vec::new();
            for mapped_base in bases {
                let mut gt = String::from("0");
                if *mapped_base != ref_base {
                    if *mapped_base == b'-' {
                        gt = ".".to_string();
                    } else {
                        let alt_base = match *mapped_base {
                            b'A' => Base::A,
                            b'C' => Base::C,
                            b'G' => Base::G,
                            b'T' => Base::T,
                            _ => Base::N,
                        };
                        if !alt_bases.contains(&alt_base) {
                            alt_bases.push(alt_base);
                        }
                        gt = (alt_bases.iter().position(|&r| r == alt_base).unwrap() + 1).to_string();
                    }
                }
                let field = Field::new(Key::Genotype, Some(Value::String(gt)));
                genotype_vec.push(Genotype::try_from(vec![field]).expect("Could not construct genotypes"));
            }
            if alt_bases.len() > 0 {
                let genotypes = Genotypes::new(
                    keys.clone(),
                    genotype_vec,
                );
                let alt_alleles: Vec<Allele> = alt_bases.iter().map(|a| Allele::Bases(vec![*a])).collect();
                let record = vcf::Record::builder()
                    .set_chromosome(self.chrom_names[*map_chrom].parse().expect("Invalid chromosome name"))
                    .set_position(Position::from(*map_pos + 1))
                    .add_reference_base(ref_allele)
                    .set_alternate_bases(AlternateBases::from(alt_alleles))
                    .set_genotypes(genotypes)
                    .build().expect("Could not construct record");
                writer.write_record(&record)?;
            }
            next_pos = *map_pos + 1;
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
        let half_split_len = (self.k - 1) / 2;
        for (sample_idx, sample_name) in self.mapped_names.iter().enumerate() {
            let sample_vars = self.mapped_variants.slice(s![.., sample_idx]);
            let mut seq: Vec<u8> = Vec::new();
            seq.reserve(self.total_size);

            // TODO: test with INDELs
            // if this proves tricky, it might be better to iterate through non-missing
            // matches and paste each flanking end and middle base into the right place
            // in the vec (and skip to next match where right end is beyond last paste)
            let (mut next_pos, mut curr_chrom) = (0, 0);
            for ((map_chrom, map_pos), base) in self.mapped_pos.iter().zip(sample_vars.iter()) {
                // Move forward to next chromosome/contig
                if *map_chrom > curr_chrom {
                    seq.extend_from_slice(&self.seq[curr_chrom][next_pos..(next_pos + half_split_len)]);
                    seq.extend(vec![b'-'; self.seq[curr_chrom].len() - (next_pos + half_split_len)]);
                    curr_chrom += 1;
                    next_pos = 0;
                }
                if *map_pos > next_pos {
                    // Missing bases, if no k-mers mapped over a region
                    seq.extend(vec![b'-'; *map_pos - next_pos - half_split_len]);
                    seq.extend_from_slice(&self.seq[curr_chrom][(*map_pos - half_split_len)..*map_pos]);
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
            seq.extend(vec![b'-'; self.seq[curr_chrom].len() - (next_pos + half_split_len)]);
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
