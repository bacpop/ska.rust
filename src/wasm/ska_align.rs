//! Library that implements sequence alignment for a WebAssembly environment

use crate::logw;
use crate::ska_dict::bit_encoding::UInt;
use crate::wasm::alignment_matrix::{
    strict_triangle_index, strict_triangle_len, AlignmentExportSession, AlignmentMatrix,
};

use speedytree::DistanceMatrix;
use speedytree::{Canonical, NeighborJoiningSolver};

#[cfg(target_family = "wasm")]
#[derive(Debug)]
/// Main struct for alignment in a WebAssembly environment.
///
/// Distances are held once in strict lower-triangle row order. An imported
/// alignment has distances and names but no sparse split-kmer matrix.
pub struct SkaAlign<IntT> {
    /// Sparse split-kmer postings, absent for an imported alignment.
    matrix: Option<AlignmentMatrix<IntT>>,
    /// Pair order is (1,0), (2,0), (2,1), …; no diagonal entries are stored.
    distances: Vec<u32>,
    /// Sample names in insertion order.
    names: Vec<String>,
    /// True if at least one imported pairwise distance was capped at u32::MAX.
    distances_capped: bool,
    /// Resumable split-kmer alignment exporter.
    export_session: Option<AlignmentExportSession>,
    export_samples_processed: usize,
    export_part_complete: bool,
    export_complete: bool,
}

impl<IntT> SkaAlign<IntT>
where
    IntT: for<'a> UInt<'a> + 'static,
    u128: num_traits::AsPrimitive<IntT>,
{
    #[cfg(target_family = "wasm")]
    /// Constructor for an alignment built from split-kmer samples.
    pub fn new() -> Self {
        Self {
            matrix: Some(AlignmentMatrix::default()),
            distances: Vec::new(),
            names: Vec::new(),
            distances_capped: false,
            export_session: None,
            export_samples_processed: 0,
            export_part_complete: false,
            export_complete: false,
        }
    }

    #[cfg(target_family = "wasm")]
    /// Constructor for an already aligned FASTA matrix. The distance vector is moved in.
    pub fn from_imported_distances(
        names: Vec<String>,
        distances: Vec<u32>,
        distances_capped: bool,
    ) -> Result<Self, String> {
        let expected = strict_triangle_len(names.len())
            .ok_or_else(|| "Imported alignment has too many samples.".to_string())?;
        if distances.len() != expected {
            return Err("Imported distance count does not match the sample count.".into());
        }
        Ok(Self {
            matrix: None,
            distances,
            names,
            distances_capped,
            export_session: None,
            export_samples_processed: 0,
            export_part_complete: false,
            export_complete: false,
        })
    }

    #[cfg(target_family = "wasm")]
    /// Stores packed split-kmer entries and updates this alignment's distance triangle in place.
    pub(crate) fn add_sample_words(
        &mut self,
        name: &str,
        words: Vec<u32>,
        on_stage: &dyn Fn(&str),
    ) -> Result<(), String> {
        on_stage("distances");
        let row_start = self.distances.len();
        let row_end = row_start
            .checked_add(self.names.len())
            .ok_or_else(|| "Distance matrix size overflowed.".to_string())?;
        if row_end > self.distances.capacity() {
            return Err("Distance capacity was not reserved for the planned batch.".into());
        }

        let matrix = self
            .matrix
            .as_mut()
            .ok_or_else(|| "Imported alignments cannot accept split-kmer samples.".to_string())?;
        self.distances.resize(row_end, 0);
        if let Err(error) = matrix.add_sample_words(words, &mut self.distances[row_start..row_end])
        {
            self.distances.truncate(row_start);
            return Err(error);
        }

        logw(
            &format!(
                "Computed distances for sample {} and stored it in the alignment",
                name
            ),
            Some("info"),
        );
        self.names.push(name.to_string());
        Ok(())
    }

    #[cfg(target_family = "wasm")]
    /// Builds a tree using precomputed pairwise distances and explicit tip labels.
    pub fn build_tree(
        &self,
        tree_tip_names: &[String],
        on_stage: &dyn Fn(&str),
    ) -> Result<String, String> {
        let n = self.names.len();
        if n < 3 || tree_tip_names.len() != n {
            return Err("Tree construction requires at least three matching sample names.".into());
        }
        let expected =
            strict_triangle_len(n).ok_or_else(|| "Distance matrix size overflowed.".to_string())?;
        if self.distances.len() != expected {
            return Err("Distance count does not match the sample count.".into());
        }

        logw(
            &format!("Initiating alignment in SkaAlign with {} input files.", n),
            None,
        );
        on_stage("tree-distances");
        logw("Building SpeedyTree distance matrix", None);

        let mut names = Vec::new();
        names
            .try_reserve_exact(n)
            .map_err(|error| format!("Could not allocate tree sample names: {error}"))?;
        for name in tree_tip_names {
            names.push(name.clone());
        }

        let mut matrix = Vec::new();
        matrix
            .try_reserve_exact(n)
            .map_err(|error| format!("Could not allocate tree matrix rows: {error}"))?;
        for row in 0..n {
            let mut values = Vec::new();
            values
                .try_reserve_exact(n)
                .map_err(|error| format!("Could not allocate tree matrix row: {error}"))?;
            for column in 0..n {
                let distance = if row == column {
                    0
                } else {
                    let (lower_row, lower_column) = if row > column {
                        (row, column)
                    } else {
                        (column, row)
                    };
                    let index = strict_triangle_index(lower_row, lower_column)
                        .ok_or_else(|| "Invalid distance matrix index.".to_string())?;
                    self.distances[index]
                };
                values.push(f64::from(distance));
            }
            matrix.push(values);
        }

        let distance_matrix = DistanceMatrix::build(matrix, names)
            .map_err(|error| format!("Could not build SpeedyTree distance matrix: {error}"))?;
        on_stage("tree");
        logw("Calculating tree", None);
        let tree = NeighborJoiningSolver::<Canonical>::default(distance_matrix)
            .solve()
            .map_err(|error| format!("Could not build neighbour-joining tree: {error}"))?;

        logw("Obtaining tree", None);
        Ok(speedytree::to_newick(&tree))
    }

    #[cfg(target_family = "wasm")]
    /// Borrow the strict lower-triangle distance vector.
    pub fn get_flat_distances(&self) -> &[u32] {
        &self.distances
    }

    #[cfg(target_family = "wasm")]
    /// Whether any imported pairwise distance was capped.
    pub fn distances_capped(&self) -> bool {
        self.distances_capped
    }

    #[cfg(target_family = "wasm")]
    /// Gets number of loaded samples.
    pub fn get_size(&self) -> usize {
        self.names.len()
    }

    /// Whether additions are closed or unavailable for an imported alignment.
    pub fn is_sealed(&self) -> bool {
        self.matrix
            .as_ref()
            .map_or(true, AlignmentMatrix::is_sealed)
    }

    /// Reserve all distance entries through the planned total in one call.
    pub fn set_expected_samples(&mut self, total: usize) -> Result<(), String> {
        let matrix = self
            .matrix
            .as_mut()
            .ok_or_else(|| "Imported alignments cannot accept new samples.".to_string())?;
        let required = strict_triangle_len(total)
            .ok_or_else(|| "Distance matrix size overflowed.".to_string())?;
        if required < self.distances.len() {
            return Err("Planned sample count is below the current alignment size.".into());
        }
        if self.distances.capacity() < required {
            self.distances
                .try_reserve_exact(required - self.distances.len())
                .map_err(|error| format!("Could not reserve distance matrix: {error}"))?;
        }
        matrix.set_expected_samples(total)?;
        Ok(())
    }

    /// Whether this handle has split-kmer sites that can be exported as FASTA.
    pub fn can_export_alignment(&self) -> bool {
        self.matrix.is_some()
    }

    /// Start a resumable export. Each returned chunk is one gzip member.
    pub fn begin_alignment_export(&mut self, samples_per_part: usize) -> Result<(), String> {
        if self.export_session.is_some() {
            return Err("An alignment export is already in progress.".into());
        }
        let matrix = self.matrix.as_mut().ok_or_else(|| {
            "Imported alignments cannot be exported as split-kmer FASTA.".to_string()
        })?;
        self.export_session = Some(matrix.begin_export(&self.names, samples_per_part)?);
        self.export_samples_processed = 0;
        self.export_part_complete = false;
        self.export_complete = false;
        Ok(())
    }

    /// Compress and return the next batch, retaining only the reusable packed workspace.
    pub fn next_alignment_export_chunk(
        &mut self,
        on_sample: &dyn Fn(usize),
    ) -> Result<Vec<u8>, String> {
        let session = self
            .export_session
            .as_mut()
            .ok_or_else(|| "No alignment export is in progress.".to_string())?;
        let matrix = self
            .matrix
            .as_mut()
            .ok_or_else(|| "No extracted alignment is available.".to_string())?;
        let chunk = matrix.next_export_chunk(&self.names, session, on_sample)?;
        self.export_samples_processed = chunk.samples_processed;
        self.export_part_complete = chunk.part_complete;
        self.export_complete = chunk.complete;
        if chunk.complete {
            self.export_session = None;
        }
        Ok(chunk.bytes)
    }

    pub fn alignment_export_samples_processed(&self) -> usize {
        self.export_samples_processed
    }

    pub fn alignment_export_part_complete(&self) -> bool {
        self.export_part_complete
    }

    pub fn alignment_export_complete(&self) -> bool {
        self.export_complete
    }

    /// Release the frozen sparse sites after the caller has accepted the final output.
    pub fn finish_alignment_export(&mut self) -> Result<(), String> {
        if !self.export_complete {
            return Err("The alignment export is not complete yet.".into());
        }
        self.matrix
            .as_mut()
            .ok_or_else(|| "No extracted alignment is available.".to_string())?
            .release_export()
    }

    /// Discard exporter state and reset cursors for a retry from sample zero.
    pub fn cancel_alignment_export(&mut self) {
        if let Some(matrix) = self.matrix.as_mut() {
            matrix.reset_export_cursors();
        }
        self.export_session = None;
        self.export_samples_processed = 0;
        self.export_part_complete = false;
        self.export_complete = false;
    }
}

/// Produce Newick-safe, unique labels while preserving the original sample order.
#[cfg(target_family = "wasm")]
pub(crate) fn make_tree_tip_names(names: &[String]) -> Vec<String> {
    use std::collections::{HashMap, HashSet};

    let stems: Vec<String> = names
        .iter()
        .map(|name| {
            let lower = name.to_ascii_lowercase();
            let end = if lower.ends_with(".gz") {
                name.len() - 3
            } else {
                name.len()
            };
            let suffixes = [
                ".fa", ".fna", ".fasta", ".fas", ".fq", ".fnq", ".fastq", ".aln",
            ];
            let stem_end = suffixes
                .iter()
                .find(|suffix| lower[..end].ends_with(**suffix))
                .map_or(name.len(), |suffix| end - suffix.len());
            let stem: String = name[..stem_end]
                .chars()
                .map(|c| {
                    if c.is_ascii_alphanumeric() || matches!(c, '.' | '_' | '-') {
                        c
                    } else {
                        '_'
                    }
                })
                .collect();
            if stem.is_empty() {
                "sample".to_string()
            } else {
                stem
            }
        })
        .collect();
    // Reserve natural stems so a duplicate "a" cannot take another sample's "a__2".
    let mut reserved: HashSet<String> = stems.iter().cloned().collect();
    let mut used = HashSet::new();
    let mut next_suffix: HashMap<String, usize> = HashMap::new();
    stems
        .into_iter()
        .map(|stem| {
            if used.insert(stem.clone()) {
                return stem;
            }
            let next = next_suffix.entry(stem.clone()).or_insert(2);
            loop {
                let candidate = format!("{stem}__{next}");
                *next += 1;
                if reserved.insert(candidate.clone()) {
                    return candidate;
                }
            }
        })
        .collect()
}
