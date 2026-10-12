//! Utilities and ports for compiling to WebAssembly targets
pub(crate) mod alignment_matrix;
#[cfg(target_family = "wasm")]
pub mod fastx_wasm;
#[cfg(target_family = "wasm")]
pub mod ska_align;
#[cfg(target_family = "wasm")]
pub mod ska_cluster;
#[cfg(target_family = "wasm")]
pub mod ska_map;
