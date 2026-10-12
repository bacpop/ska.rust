//! Clustering of samples by SNP distance using a graph-based approach.
//!
//! Builds an undirected graph where edges connect samples within a SNP distance
//! threshold, then extracts connected components as clusters.

use std::cmp::Reverse;
#[cfg(not(target_family = "wasm"))]
use std::io::Write;

use hashbrown::HashMap;
#[cfg(not(target_family = "wasm"))]
use petgraph::dot::{Config, Dot};
use petgraph::graph::NodeIndex;
use petgraph::visit::Bfs;
use petgraph::Graph;

#[cfg(any(target_family = "wasm", test))]
use crate::wasm::alignment_matrix::{strict_triangle_index, strict_triangle_len};
use crate::merge_ska_array::VariantDist;

type ClusterResults = (
    HashMap<String, usize>,
    Graph<String, (), petgraph::Undirected>,
);

/// Shared graph-building and component-labelling logic.
fn build_clusters(names: &[String], connected: impl Fn(usize, usize) -> bool) -> ClusterResults {
    let mut graph = Graph::new_undirected();
    let n = names.len();
    for name in names {
        graph.add_node(name.clone());
    }
    for i in 0..n {
        for j in (i + 1)..n {
            if connected(i, j) {
                graph.add_edge(NodeIndex::new(i), NodeIndex::new(j), ());
            }
        }
    }

    // Label connected components via BFS traversal
    let mut component = vec![usize::MAX; names.len()];
    let mut comp_id = 0usize;
    for start in graph.node_indices() {
        if component[start.index()] == usize::MAX {
            let mut bfs = Bfs::new(&graph, start);
            while let Some(node) = bfs.next(&graph) {
                component[node.index()] = comp_id;
            }
            comp_id += 1;
        }
    }

    // Sort components by size descending, assign 1-indexed IDs
    let mut comp_sizes: Vec<(usize, usize)> = (0..comp_id)
        .map(|c| (c, component.iter().filter(|&&x| x == c).count()))
        .collect();
    comp_sizes.sort_unstable_by_key(|entry| Reverse(entry.1));
    let mut id_map = vec![0usize; comp_id];
    for (rank, (orig_id, _)) in comp_sizes.iter().enumerate() {
        id_map[*orig_id] = rank + 1;
    }

    // Build result map
    let mut result: HashMap<String, usize> = HashMap::new();
    for (i, name) in names.iter().enumerate() {
        result.insert(name.clone(), id_map[component[i]]);
    }

    (result, graph)
}

/// Build clusters from pairwise distances.
///
/// Returns a map of sample name → 1-indexed cluster ID (largest cluster = 1)
/// and the graph used for clustering.
pub fn cluster_distances(
    names: &[String],
    distances: &[Vec<VariantDist>],
    threshold: f64,
) -> ClusterResults {
    build_clusters(names, |i, j| {
        distances[i][j - i - 1].distance() <= threshold
    })
}

/// Convert the public floating-point threshold to the integer count domain.
/// Fractional values retain the existing `distance <= threshold` semantics.
#[cfg(any(target_family = "wasm", test))]
pub(crate) fn checked_distance_cutoff(threshold: f64) -> Result<u32, String> {
    if !threshold.is_finite() || threshold < 0.0 || threshold > f64::from(u32::MAX) {
        return Err("SNP threshold must be finite and between 0 and u32::MAX.".into());
    }
    Ok(threshold.floor() as u32)
}

/// Build clusters from strict lower-triangle integer distances (WASM entry point).
#[cfg(any(target_family = "wasm", test))]
pub fn cluster_distances_flat(
    names: &[String],
    flat: &[u32],
    threshold: u32,
) -> Result<ClusterResults, String> {
    let expected = strict_triangle_len(names.len())
        .ok_or_else(|| "Sample count exceeds the distance index range.".to_string())?;
    if flat.len() != expected {
        return Err("Distance count does not match the sample count.".into());
    }
    Ok(build_clusters(names, |i, j| {
        let index =
            strict_triangle_index(j, i).expect("cluster pairs are strictly above the diagonal");
        flat[index] <= threshold
    }))
}

/// Write cluster CSV and DOT graph files.
#[cfg(not(target_family = "wasm"))]
pub fn write_graph(
    graph: &Graph<String, (), petgraph::Undirected>,
    cluster_map: &HashMap<String, usize>,
    output_prefix: &Option<String>,
) {
    let prefix = output_prefix.as_deref().unwrap_or("ska_dist_clusters");

    // Write clusters CSV
    let csv_path = format!("{prefix}.clusters.csv");
    let mut csv = std::fs::File::create(&csv_path).expect("Cannot create clusters CSV");
    writeln!(csv, "id,Cluster__autocolour").unwrap();
    let mut entries: Vec<(&String, &usize)> = cluster_map.iter().collect();
    entries.sort_by_key(|(name, &id)| (id, name.as_str().to_owned()));
    for (name, id) in entries {
        writeln!(csv, "{name},{id}").unwrap();
    }
    log::info!("Written clusters to {csv_path}");

    // Write DOT file
    let dot_path = format!("{prefix}.graph.dot");
    let dot = Dot::with_config(graph, &[Config::EdgeNoLabel]);
    let mut dot_file = std::fs::File::create(&dot_path).expect("Cannot create DOT file");
    write!(dot_file, "{:?}", dot).unwrap();
    log::info!("Written network to {dot_path}");
}

#[cfg(test)]
mod tests {
    use super::{checked_distance_cutoff, cluster_distances_flat};

    #[test]
    fn strict_lower_triangle_clustering_uses_every_pair_in_order() {
        let names = ["a", "b", "c", "d"].map(str::to_owned);
        // Pair values: a-b=1, a-c=2, b-c=3, a-d=4, b-d=5, c-d=6.
        let distances = [1, 2, 3, 4, 5, 6];
        let (clusters, graph) = cluster_distances_flat(&names, &distances, 2).unwrap();

        assert_eq!(graph.edge_count(), 2);
        assert_eq!(clusters["a"], clusters["b"]);
        assert_eq!(clusters["a"], clusters["c"]);
        assert_ne!(clusters["a"], clusters["d"]);
    }

    #[test]
    fn public_threshold_conversion_preserves_fractional_cutoff_semantics() {
        assert_eq!(checked_distance_cutoff(2.0).unwrap(), 2);
        assert_eq!(checked_distance_cutoff(2.9).unwrap(), 2);
        assert!(checked_distance_cutoff(-0.1).is_err());
        assert!(checked_distance_cutoff(f64::NAN).is_err());
        assert!(checked_distance_cutoff(f64::INFINITY).is_err());
        assert!(checked_distance_cutoff(f64::from(u32::MAX) + 1.0).is_err());
    }
}
