//! Library that implements sequence alignment for a WebAssembly environment

use hashbrown::HashMap;

use crate::logw;
use crate::ska_dict::bit_encoding::UInt;
use crate::ska_dict::SkaDict;
use crate::QualFilter;
use crate::QualOpts;

use speedytree::DistanceMatrix;
use speedytree::{Canonical, NeighborJoiningSolver};

#[cfg(target_family = "wasm")]
#[derive(Debug, Clone, Default)]
/// Main struct for alignment in a WebAssembly environment.
///
/// Distances are computed incrementally as files are added with [`add_file()`],
/// so only the raw split-kmer HashMaps (one per past sample) stay in memory —
/// not the full [`SkaDict`] objects with their bloom filter buffers.
pub struct SkaAlign<IntT> {
    /// Raw split-kmer maps for all samples loaded so far
    past_kmers: Vec<HashMap<IntT, u8>>,
    /// Pairwise distance matrix in lower-triangular row storage: `distances[i]`
    /// has `i + 1` entries and `distances[i][j]` (j ≤ i) = SNP distance between
    /// sample i and sample j. Use `get_flat_distances()` for upper-triangle order.
    distances: Vec<Vec<u32>>,
    /// Sample names in insertion order
    names: Vec<String>,
    /// k-value being used
    k: usize,
    /// Use canonical k-mers or not
    rc: bool,
}

impl<IntT> SkaAlign<IntT>
where
    IntT: for<'a> UInt<'a>,
{
    #[cfg(target_family = "wasm")]
    /// Constructor of a SkaAlign struct
    pub fn new(k: usize, rc: bool) -> Self {
        Self {
            past_kmers: Vec::new(),
            distances: Vec::new(),
            names: Vec::new(),
            k,
            rc,
        }
    }

    #[cfg(target_family = "wasm")]
    /// Adds a file, builds its split-kmer dict, computes distances against all
    /// previously loaded samples, then retains only the raw kmer map.
    ///
    /// This keeps peak memory at O(N × dict_size) rather than O(N × SkaDict_size),
    /// saving the ~24 MB bloom filter buffer per FASTQ sample.
    pub fn add_file(
        &mut self,
        file1: &web_sys::File,
        file2: Option<&web_sys::File>,
        proportions_reads: Option<f64>,
        min_count: u16,
        min_qual: u8,
        qual_filter: QualFilter,
        name: &str,
        idx: usize,
    ) {
        let ska_dict = SkaDict::new(
            self.k,
            idx,
            (file1, file2),
            name,
            self.rc,
            &QualOpts {
                min_count,
                min_qual,
                qual_filter,
            },
            proportions_reads,
        );

        // Compute distances against all already-loaded samples (lower triangle)
        let n_prev = self.past_kmers.len();
        let mut row = vec![0u32; n_prev + 1]; // last entry = self-distance (0)
        for (j, prev_kmers) in self.past_kmers.iter().enumerate() {
            let mut dist = 0u32;
            for (kmer, base) in ska_dict.kmers().iter() {
                if let Some(prev_base) = prev_kmers.get(kmer) {
                    if prev_base != base {
                        dist += 1;
                    }
                }
            }
            row[j] = dist;
        }
        self.distances.push(row);
        self.names.push(name.to_string());

        // Keep only the kmer map — drops the bloom buffer and SkaDict wrapper
        self.past_kmers.push(ska_dict.into_kmers());
    }

    #[cfg(target_family = "wasm")]
    /// Performs the alignment using precomputed pairwise distances.
    pub fn align(&mut self, file_names: &[String]) -> String {
        let n = self.names.len();
        logw(
            &format!("Initiating alignment in SkaAlign with {} input files.", n,),
            None,
        );

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

        logw(&format!("{:?}", phylip_format), None);
        logw("Converting matrix to DistanceMatrix struct.", None);

        let d = DistanceMatrix::read_from_phylip(phylip_format.as_bytes()).unwrap();

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

    #[cfg(target_family = "wasm")]
    /// Iterate over (sample_index, sample_name, kmer_map) for all loaded samples.
    ///
    /// Used to build a [`crate::merge_ska_dict::MergeSkaDict`] for FASTA output
    /// without needing the full [`crate::ska_dict::SkaDict`] objects.
    pub fn iter_kmers(&self) -> impl Iterator<Item = (usize, &str, &HashMap<IntT, u8>)> {
        self.past_kmers
            .iter()
            .enumerate()
            .map(|(i, km)| (i, self.names[i].as_str(), km))
    }
}
