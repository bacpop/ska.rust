//! Library that implements sequence alignment for a WebAssembly environment

use crate::alignment_matrix::AlignmentMatrix;
use crate::logw;
use crate::ska_dict::bit_encoding::UInt;

use speedytree::DistanceMatrix;
use speedytree::{Canonical, NeighborJoiningSolver};

#[cfg(target_family = "wasm")]
#[derive(Debug, Clone, Default)]
/// Main struct for alignment in a WebAssembly environment.
///
/// Distances are computed incrementally as packed sample entries are added.
/// The matrix retains one shared k-mer-to-row map and sample-ID/base packed calls,
/// rather than one raw split-kmer map per sample.
pub struct SkaAlign<IntT> {
    /// Sparse split-kmer postings.
    matrix: AlignmentMatrix<IntT>,
    /// Pairwise distance matrix in lower-triangular row storage: `distances[i]`
    /// has `i + 1` entries and `distances[i][j]` (j ≤ i) = SNP distance between
    /// sample i and sample j. Use `get_flat_distances()` for upper-triangle order.
    distances: Vec<Vec<u32>>,
    /// Sample names in insertion order
    names: Vec<String>,
}

impl<IntT> SkaAlign<IntT>
where
    IntT: for<'a> UInt<'a>,
{
    #[cfg(target_family = "wasm")]
    /// Constructor of a SkaAlign struct
    pub fn new() -> Self {
        Self {
            matrix: AlignmentMatrix::default(),
            distances: Vec::new(),
            names: Vec::new(),
        }
    }

    #[cfg(target_family = "wasm")]
    /// Stores packed split-kmer entries and computes their distance row in one pass.
    pub fn add_sample_words(
        &mut self,
        name: &str,
        words: Vec<u32>,
        words_per_entry: usize,
        on_stage: &dyn Fn(&str),
    ) -> Result<(), String> {
        on_stage("distances");
        let row = self.matrix.add_sample_words(words, words_per_entry)?;
        logw(
            &format!(
                "Computed distances for sample {} and stored it in the sparse alignment matrix",
                name
            ),
            Some("info"),
        );

        self.distances.push(row);
        self.names.push(name.to_string());
        Ok(())
    }

    #[cfg(target_family = "wasm")]
    /// Performs the alignment using precomputed pairwise distances.
    pub fn align(&mut self, file_names: &[String], on_stage: &dyn Fn(&str)) -> String {
        let n = self.names.len();
        logw(
            &format!("Initiating alignment in SkaAlign with {} input files.", n,),
            None,
        );

        on_stage("tree-distances");
        logw("Creating pairwise distances matrix as text.", None);

        let mut phylip_format = format!("{}\n", n);

        for i in 0..n {
            phylip_format += file_names[i]
                .to_string()
                .replace(" ", "_")
                .replace(".fasta", "")
                .replace(".fa", "")
                .replace(".fastq", "")
                .replace(".fq", "")
                .as_str();
            for j in 0..n {
                // distances matrix is lower-triangular: distances[i][j] exists for j <= i
                let dist = if i == j {
                    0u32
                } else if j < i {
                    self.distances[i][j]
                } else {
                    self.distances[j][i]
                };
                phylip_format += &format!("\t{dist}");
            }
            phylip_format += "\n";
        }

        logw(
            &format!("Prepared distance matrix for {} samples", n),
            Some("info"),
        );

        let d = DistanceMatrix::read_from_phylip(phylip_format.as_bytes()).unwrap();

        on_stage("tree");
        logw("Calculating tree", None);
        let tree = NeighborJoiningSolver::<Canonical>::default(d.clone())
            .solve()
            .unwrap();

        logw("Obtaining tree", None);

        speedytree::to_newick(&tree)
    }

    #[cfg(target_family = "wasm")]
    /// Returns pairwise distances as a flat `Vec<f64>` in upper-triangle row-major
    /// order: all pairs (i, j) with i < j, i increasing from 0. The distance
    /// between sample i and sample j (i < j) is at index `i*n - i*(i+1)/2 + (j-i-1)`.
    /// Internally reads from the lower-triangular store as `distances[j][i]`.
    pub fn get_flat_distances(&self) -> Vec<f64> {
        let n = self.names.len();
        let mut flat = Vec::with_capacity(n * (n - 1) / 2);
        for i in 0..n {
            for j in (i + 1)..n {
                flat.push(self.distances[j][i] as f64);
            }
        }
        flat
    }

    #[cfg(target_family = "wasm")]
    /// Gets number of loaded samples
    pub fn get_size(&self) -> usize {
        self.names.len()
    }

    /// Whether alignment export has permanently closed sample additions.
    pub fn is_sealed(&self) -> bool {
        self.matrix.is_sealed()
    }

    /// Close additions and lazily stream the final sparse array to gzip.
    pub fn export_fasta_gzip(&mut self, on_sample: &dyn Fn(usize)) -> Result<Vec<u8>, String> {
        self.matrix.export_fasta_gzip(&self.names, on_sample)
    }
}
