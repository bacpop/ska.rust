//! Clustering of samples by SNP distance using a graph-based approach.
//!
//! Builds an undirected graph where edges connect samples within a SNP distance
//! threshold, then extracts connected components as clusters.

#[cfg(not(target_family = "wasm"))]
use std::io::Write;

use hashbrown::HashMap;
#[cfg(not(target_family = "wasm"))]
use petgraph::dot::{Config, Dot};
use petgraph::graph::NodeIndex;
use petgraph::visit::Bfs;
use petgraph::Graph;

use crate::merge_ska_array::VariantDist;

/// Shared graph-building and component-labelling logic.
fn build_clusters(
    names: &[String],
    threshold: f64,
    dist_fn: impl Fn(usize, usize) -> f64,
) -> (
    HashMap<String, usize>,
    Graph<String, (), petgraph::Undirected>,
) {
    let mut graph = Graph::new_undirected();
    let n = names.len();
    for name in names {
        graph.add_node(name.clone());
    }
    for i in 0..n {
        for j in (i + 1)..n {
            if dist_fn(i, j) <= threshold {
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
    comp_sizes.sort_unstable_by(|a, b| b.1.cmp(&a.1));
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
) -> (
    HashMap<String, usize>,
    Graph<String, (), petgraph::Undirected>,
) {
    build_clusters(names, threshold, |i, j| distances[i][j - i - 1].distance())
}

/// Build clusters from a flat upper-triangle distance array (wasm entry point).
#[cfg(target_family = "wasm")]
pub fn cluster_distances_flat(
    names: &[String],
    flat: &[f64],
    threshold: f64,
) -> (
    HashMap<String, usize>,
    Graph<String, (), petgraph::Undirected>,
) {
    let n = names.len();
    build_clusters(names, threshold, |i, j| {
        let idx = i * n - i * (i + 1) / 2 + (j - i - 1);
        flat[idx]
    })
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
