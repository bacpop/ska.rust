//! Library that implements mapping sequences for a WebAssembly environment

use crate::ska_dict::SkaDict;
use crate::ska_dict::bit_encoding::{UInt, RC_IUPAC};
use crate::ska_ref::RefSka;
use crate::QualFilter;
use crate::QualOpts;
use std::marker::PhantomData;

#[cfg(feature = "wasm")]
#[derive(Debug, Clone, Default)]
/// Variant structure
pub struct Variant {
    /// Chromosome where the variant is
    pub chrom: usize,
    /// Position in the sequence
    pub pos: usize,
    /// Base of the variant
    pub base: u8,
}

#[cfg(feature = "wasm")]
#[derive(Debug, Clone, Default)]
/// SkaMap structure, to hold both mapping information as well as creation/usage methods
pub struct SkaMap<IntT> {
    /// Positions of mapped bases as (chrom, pos, base)
    mapped_bases: Vec<Variant>,
    /// Phantom attribute for having SkaMap variations depending on generics
    phantom: PhantomData<IntT>,
}

impl<IntT> SkaMap<IntT>
where
    IntT: for<'a> UInt<'a>,
{
    #[cfg(feature = "wasm")]
    /// Constructor for the SkaMap struct.
    pub fn new(
        reference: &RefSka<IntT>,
        file1: &web_sys::File,
        file2: Option<&web_sys::File>,
        proportions_reads: Option<f64>,
        rc: bool,
    ) -> Self {
        // TODO - two files. Should just be able to add an 'add more k-mers' method on the struct to accept second file if given

        let qualities = QualOpts {
            min_count: 1,
            min_qual: 0,
            qual_filter: QualFilter::NoFilter,
        };

        let query_ska = SkaDict::<IntT>::new(
            reference.get_k(),
            0,
            (file1, file2),
            "",
            rc,
            &qualities,
            proportions_reads,
        );

        let mut mapped_bases = Vec::new();
        for ref_kmer in reference.kmers() {
            if let Some(kmer_match) = query_ska.kmers().get(&ref_kmer.kmer) {
                if ref_kmer.rc {
                    mapped_bases.push(Variant {
                        chrom: ref_kmer.chrom,
                        pos: ref_kmer.pos,
                        base: RC_IUPAC[*kmer_match as usize],
                    })
                } else {
                    mapped_bases.push(Variant {
                        chrom: ref_kmer.chrom,
                        pos: ref_kmer.pos,
                        base: *kmer_match,
                    })
                }
            }
        }
        Self {
            mapped_bases,
            phantom: PhantomData,
        }
    }

    #[cfg(feature = "wasm")]
    /// Gets the mapped bases attribute.
    pub fn get_mapped_bases(&self) -> &Vec<Variant> {
        &self.mapped_bases
    }
}
