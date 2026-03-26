//! Clustering adapter for the WebAssembly environment

use wasm_bindgen::prelude::*;

use crate::cluster::cluster_distances_flat;
use crate::AlignData;

#[wasm_bindgen]
/// Cluster samples by SNP distance threshold using a pre-computed `AlignData`.
///
/// Call this after `AlignData.align()` has been run.
/// `threshold` – SNP distance cutoff
/// Returns JSON string: `{"sample_name": cluster_id, ...}` (1-indexed, largest cluster = 1).
pub fn ska_cluster(data: &AlignData, threshold: f64) -> String {
    let names = data.names();
    let flat = data.flat_distances();
    let (cluster_map, _graph) = cluster_distances_flat(names, flat, threshold);
    let mut result = json::JsonValue::new_object();
    for (name, id) in &cluster_map {
        result[name] = (*id).into();
    }
    result.dump()
}
