
use hashbrown::HashMap;
extern crate needletail;
use ndarray::Array1;
use needletail::{parse_fastx_file, parser::Format};
use libm::lgamma;
use cached::proc_macro::cached;

use crate::ska_dict::bit_encoding::UInt;

pub struct CoverageHistogram<IntT> {
    /// Dictionary of k-mers
    kmer_counts: HashMap<IntT, u32>,
    /// Count histogram
    counts: Vec<u32>,
    /// K-mer size
    k: usize,
    /// Whether reverse complement split k-mers were used
    rc: bool,
}

impl<IntT> CoverageHistogram<IntT>
where
    IntT: for<'a> UInt<'a>,
{
    // e.g. f = |x: &Array1<f64>| log_likelihood(x, self.counts);

}

// Called by lib.rs
pub fn coverage_cutoff<IntT: for<'a> UInt<'a>>(fastq1: &String, fastq2: &String, k: usize, rc: bool, threads: usize) {

}

// log-sum-exp
fn lse(a: f64, b: f64) -> f64 {
    let xstar = f64::max(a, b);
    xstar + f64::ln(f64::exp(a - xstar) + f64::exp(b - xstar))
}

// TODO: 64 needs to be cached as bytes

// Natural log of Poisson density
#[cached]
fn ln_dpois(x: u64, lambda: [u8; 8]) -> f64 {
    let lambda_f = f64::from_le_bytes(lambda);
    x as f64 * (lambda_f) - lgamma(x as f64 + 1.0) - lambda_f
}

// error component
#[cached]
fn a(w0: f64, i: f64) -> f64 {
    f64::ln(w0) + ln_dpois(i, 1.0)
}

// coverage component
#[cached]
fn b(w0: f64, c: f64, i: f64) -> f64 {
    f64::ln(1.0 - w0) + ln_dpois(i, c)
}

// Mixture likelihood
fn log_likelihood(pars: &Array1<f64>, counts: &[f64]) -> f64 {
    let w0 = pars[0];
    let c = pars[1];
    let mut ll = 0.0;
    if w0 > 1.0 || w0 < 0.0 || c < 1.0 {
        ll = f64::NEG_INFINITY;
    } else {
        for (i, c) in counts.iter().enumerate() {
            ll += counts[i] * lse(a(w0, i as f64), b(w0, *c, i as f64));
        }
    }
    ll
}

fn grad_ll(pars: &Array1<f64>, counts: &[f64]) -> Array1<f64> {
    let w0 = pars[0];
    let c = pars[1];

    let mut grad_w0 = 0.0;
    let mut grad_c = 0.0;
    for (i, c) in counts.iter().enumerate() {
        let i_f64 = i as f64;
        let a_val = a(w0, i_f64);
        let b_val = b(w0, *c, i_f64);
        let dlda = 1.0 / (1.0 + f64::exp(b_val - a_val));
        let dldb = 1.0 / (1.0 + f64::exp(a_val - b_val));
        grad_w0 += counts[i] * (dlda/w0 - dldb/(1.0 - w0));
        grad_c += counts[i] * (dldb*(i_f64/c - 1.0));
    }
    Array1::from_vec(vec![grad_w0, grad_c])
}
