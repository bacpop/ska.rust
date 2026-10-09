//! Clustering adapter for the WebAssembly environment

use wasm_bindgen::prelude::*;

use crate::cluster::{checked_distance_cutoff, cluster_distances_flat};
use crate::AlignData;

#[wasm_bindgen]
/// Cluster samples by SNP distance threshold using a pre-computed `AlignData`.
///
/// Call after `AlignData.finish_alignment()` or after importing aligned FASTA.
/// `threshold` – SNP distance cutoff
/// Returns JSON string: `{"sample_name": cluster_id, ...}` (1-indexed, largest cluster = 1).
pub fn ska_cluster(data: &AlignData, threshold: f64) -> Result<String, JsValue> {
    let cutoff = checked_distance_cutoff(threshold).map_err(|error| JsValue::from_str(&error))?;
    let names = data.names();
    let flat = data.flat_distances()?;
    let (cluster_map, _graph) =
        cluster_distances_flat(names, flat, cutoff).map_err(|error| JsValue::from_str(&error))?;
    let mut result = json::JsonValue::new_object();
    for (name, id) in &cluster_map {
        result[name] = (*id).into();
    }
    Ok(result.dump())
}
