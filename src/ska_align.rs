//! Library that implements sequence alignment for a WebAssembly environment

use crate::ska_dict::SkaDict;
use crate::ska_dict::bit_encoding::UInt;
use crate::QualFilter;
use crate::QualOpts;
use crate::logw;

use speedytree::DistanceMatrix;
use speedytree::{NeighborJoiningSolver, Canonical};

#[cfg(target_arch = "wasm32")]
#[derive(Debug, Clone, Default)]
/// Main struct for alignment in a WebAssembly environment
pub struct SkaAlign<IntT> {
    /// Positions of mapped bases as (chrom, pos)
    queries_ska: Vec<SkaDict<IntT>>,
    /// k-value being used
    k: usize,
}

impl<IntT> SkaAlign<IntT>
where
    IntT: for<'a> UInt<'a>,
{
    #[cfg(target_arch = "wasm32")]
    /// Constructor of a SkaAlign struct
    pub fn new(k: usize) -> Self {
        Self { queries_ska: Vec::new(), k }
    }

    #[cfg(target_arch = "wasm32")]
    /// Adds a file through a SkaDict.
    pub fn add_file(
        &mut self,
        file1: & web_sys::File,
        file2: Option<& web_sys::File>,
        proportions_reads: Option<f64>,
    ) {
        self.queries_ska.push(SkaDict::new(
            self.k,
            0,
            (file1, file2),
            "",
            if file2.is_some() {true} else {false},
            &QualOpts {
                min_count: 1,
                min_qual: 0,
                qual_filter: QualFilter::NoFilter,
            },
            proportions_reads
        ));
    }

    #[cfg(target_arch = "wasm32")]
    /// Performs the alignment
    pub fn align(&mut self, file_names: &Vec<String>) -> String {
        logw(&format!(
            "Initiating alignment in SkaAlign with {} input files.",
            file_names.len(),
        ), None);

        logw(&format!("Creating pairwise distances matrix as text."), None);
        let mut pairwise_distances = vec![vec![0; self.queries_ska.len()]; self.queries_ska.len()];

        let mut phylip_format = "".to_string();
        phylip_format += format!("{}\n", self.queries_ska.len()).as_str();

        for i in 0..self.queries_ska.len() {
            phylip_format += format!("{}", file_names[i]).replace(" ", "_").replace(".fasta", "").replace(".fa", "").replace(".fastq", "").replace(".fq", "").as_str();
            for j in 0..self.queries_ska.len() { // Do it on only half of the matrix
                for ref_kmer in self.queries_ska[i].kmers().iter() {
                    if let Some(kmer_match) = self.queries_ska[j].kmers().get(ref_kmer.0) {
                        if *kmer_match != *ref_kmer.1 {
                            pairwise_distances[i][j] += 1;
                        }
                    }
                }
                phylip_format += format!("\t{}", pairwise_distances[i][j]).as_str();
            }
            phylip_format += "\n";
        }

        logw(&format!("{:?}", phylip_format), None);
        logw(&format!("Converting matrix to DistanceMatrix struct."), None);

        let d = DistanceMatrix::read_from_phylip(phylip_format.as_bytes()).unwrap();

        logw(&format!("Calculating tree"), None);
        let tree = NeighborJoiningSolver::<Canonical>::default(d.clone())
            .solve()
            .unwrap();

        logw(&format!("Obtaining tree"), None);
        let newick = speedytree::to_newick(&tree);
        newick
    }

    #[cfg(target_arch = "wasm32")]
    /// Gets number of queries
    pub fn get_size(&self) -> usize {
        self.queries_ska.len()
    }
}
