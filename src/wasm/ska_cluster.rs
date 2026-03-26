//! Clustering adapter for the WebAssembly environment

use wasm_bindgen::prelude::*;

use crate::cluster::cluster_distances_flat;

#[wasm_bindgen]
/// Cluster samples by SNP distance threshold.
/// `names` – JS Array of strings
/// `distances` – flat Float64Array, upper triangle (n*(n-1)/2 entries, row-major)
/// `threshold` – SNP distance cutoff
/// Returns JSON string: {"sample_name": cluster_id, ...}
pub fn ska_cluster(
    names: js_sys::Array,
    distances: js_sys::Float64Array,
    threshold: f64,
) -> String {
    let names_vec: Vec<String> = names
        .iter()
        .map(|v| v.as_string().unwrap_or_default())
        .collect();
    let flat: Vec<f64> = distances.to_vec();
    let (cluster_map, _graph) = cluster_distances_flat(&names_vec, &flat, threshold);
    let mut result = json::JsonValue::new_object();
    for (name, id) in &cluster_map {
        result[name] = (*id).into();
    }
    result.dump()
}
