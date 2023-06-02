use core::panic;

use argmin::core::observers::{ObserverMode, SlogLogger};
use argmin::core::{
    CostFunction, Error, Executor, Gradient, State, TerminationReason::SolverConverged,
};
use argmin::solver::linesearch::condition::ArmijoCondition;
use argmin::solver::linesearch::BacktrackingLineSearch;
use argmin::solver::quasinewton::BFGS;

use libm::lgamma;

use hashbrown::HashMap;
extern crate needletail;
use needletail::{parse_fastx_file, parser::Format};

use crate::ska_dict::bit_encoding::UInt;
use crate::ska_dict::split_kmer::SplitKmer;
use crate::QualFilter;

const MAX_COUNT: usize = 1000;
const MIN_FREQ: u32 = 50;

pub struct CoverageHistogram<IntT> {
    /// K-mer size
    k: usize,
    /// Whether reverse complement split k-mers were used
    rc: bool,
    /// Dictionary of k-mers and their counts
    kmer_dict: HashMap<IntT, u32>,
    /// Count histogram
    counts: Vec<u32>,
    /// Estimated error weight
    w0: f64,
    /// Estimated coverage
    c: f64,
    /// Coverage cutoff
    cutoff: u32,
    /// Show logging
    verbose: bool,
}

impl<IntT> CoverageHistogram<IntT>
where
    IntT: for<'a> UInt<'a>,
{
    // Called by lib.rs
    pub fn new(fastq1: &String, fastq2: &String, k: usize, rc: bool, verbose: bool) -> Self {
        if !(5..=63).contains(&k) || k % 2 == 0 {
            panic!("Invalid k-mer length");
        }

        let mut cov_counts = Self {
            k,
            rc,
            kmer_dict: HashMap::default(),
            counts: vec![0; MAX_COUNT],
            w0: 0.8,
            c: 20.0,
            cutoff: 0,
            verbose,
        };

        // Check if we're working with reads first
        for fastx_file in [fastq1, fastq2] {
            let mut reader_peek = parse_fastx_file(fastx_file)
                .unwrap_or_else(|_| panic!("Invalid path/file: {}", fastx_file));
            let seq_peek = reader_peek
                .next()
                .expect("Invalid FASTA/Q record")
                .expect("Invalid FASTA/Q record");
            if seq_peek.format() != Format::Fastq {
                panic!("{fastx_file} appears to be FASTA.\nCoverage can only be used with FASTQ files, not FASTA.");
            }
        }

        log::info!("Counting k-mers");
        for fastx_file in [fastq1, fastq2] {
            let mut reader = parse_fastx_file(fastx_file)
                .unwrap_or_else(|_| panic!("Invalid path/file: {fastx_file}"));
            while let Some(record) = reader.next() {
                let seqrec = record.expect("Invalid FASTA/Q record");
                let kmer_opt = SplitKmer::new(
                    seqrec.seq(),
                    seqrec.num_bases(),
                    seqrec.qual(),
                    cov_counts.k,
                    cov_counts.rc,
                    0,
                    QualFilter::NoFilter,
                    false,
                );
                if let Some(mut kmer_it) = kmer_opt {
                    let (kmer, _base, _rc) = kmer_it.get_curr_kmer();
                    cov_counts
                        .kmer_dict
                        .entry(kmer)
                        .and_modify(|count| *count += 1)
                        .or_insert(1);
                    while let Some((kmer, _base, _rc)) = kmer_it.get_next_kmer() {
                        cov_counts
                            .kmer_dict
                            .entry(kmer)
                            .and_modify(|count| *count += 1)
                            .or_insert(1);
                    }
                }
            }
        }

        cov_counts
    }

    pub fn fit_histogram(&mut self) -> Result<u32, Error> {
        // Calculate k-mer histogram
        log::info!("Calculating k-mer histogram");
        for kmer_count in self.kmer_dict.values() {
            let kc = (*kmer_count - 1) as usize;
            if kc < MAX_COUNT {
                self.counts[kc] += 1;
            }
        }

        // Truncate count vec and covert to float
        let mut counts_f64: Vec<f64> = Vec::new();
        for hist_bin in &self.counts {
            if *hist_bin < MIN_FREQ {
                break;
            } else {
                counts_f64.push(*hist_bin as f64);
            }
        }

        log::info!("Fitting Poisson mixture model using maximum likelihood");
        let mixture_fit = MixPoisson { counts: counts_f64 };
        let init_param: Vec<f64> = vec![self.w0, self.c];
        let init_hessian: Vec<Vec<f64>> = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let linesearch = BacktrackingLineSearch::new(ArmijoCondition::new(0.0001f64)?);
        let solver = BFGS::new(linesearch);
        let mut exec = Executor::new(mixture_fit, solver).configure(|state| {
            state
                .param(init_param)
                .inv_hessian(init_hessian)
                .max_iters(100)
        });
        if self.verbose {
            exec = exec.add_observer(SlogLogger::term(), ObserverMode::Always);
        }
        let res = exec.run()?;

        // print diagnostics
        log::info!("{res}");
        if let Some(termination_reason) = res.state().get_termination_reason() {
            if *termination_reason == SolverConverged {
                // Best parameter vector
                let best = res.state().get_best_param().unwrap();
                self.w0 = best[0];
                self.c = best[1];

                // TODO calculate the coverage cutoff
                Ok(self.cutoff)
            } else {
                Err(Error::msg(format!(
                    "Optimiser did not converge: {}",
                    termination_reason.text()
                )))
            }
        } else {
            Err(Error::msg("Optimiser did not finish running"))
        }
    }

    pub fn plot_hist() {
        todo!()
    }
}

struct MixPoisson {
    counts: Vec<f64>,
}

impl CostFunction for MixPoisson {
    /// Type of the parameter vector
    type Param = Vec<f64>;
    /// Type of the return value computed by the cost function
    type Output = f64;

    /// Apply the cost function to a parameter `p`
    fn cost(&self, p: &Self::Param) -> Result<Self::Output, Error> {
        Ok(-log_likelihood(p, &self.counts))
    }
}

impl Gradient for MixPoisson {
    /// Type of the parameter vector
    type Param = Vec<f64>;
    /// Type of the gradient
    type Gradient = Vec<f64>;

    /// Compute the gradient at parameter `p`.
    fn gradient(&self, p: &Self::Param) -> Result<Self::Gradient, Error> {
        // Compute gradient of 2D Rosenbrock function
        Ok(grad_ll(p, &self.counts).iter().map(|x| -*x).collect())
    }
}

// log-sum-exp
fn lse(a: f64, b: f64) -> f64 {
    let xstar = f64::max(a, b);
    xstar + f64::ln(f64::exp(a - xstar) + f64::exp(b - xstar))
}

// Natural log of Poisson density
fn ln_dpois(x: f64, lambda: f64) -> f64 {
    x * f64::ln(lambda) - lgamma(x + 1.0) - lambda
}

// error component
fn a(w0: f64, i: f64) -> f64 {
    f64::ln(w0) + ln_dpois(i, 1.0)
}

// coverage component
fn b(w0: f64, c: f64, i: f64) -> f64 {
    f64::ln(1.0 - w0) + ln_dpois(i, c)
}

// Mixture likelihood
fn log_likelihood(pars: &[f64], counts: &[f64]) -> f64 {
    let w0 = pars[0];
    let c = pars[1];
    let mut ll = 0.0;
    if w0 > 1.0 || w0 < 0.0 || c < 1.0 {
        ll = f64::MIN;
    } else {
        for (i, count) in counts.iter().enumerate() {
            ll += *count * lse(a(w0, i as f64 + 1.0), b(w0, c, i as f64 + 1.0));
        }
    }
    ll
}

fn grad_ll(pars: &[f64], counts: &[f64]) -> Vec<f64> {
    let w0 = pars[0];
    let c = pars[1];

    let mut grad_w0 = 0.0;
    let mut grad_c = 0.0;
    for (i, count) in counts.iter().enumerate() {
        let i_f64 = i as f64 + 1.0;
        let a_val = a(w0, i_f64);
        let b_val = b(w0, c, i_f64);
        let dlda = 1.0 / (1.0 + f64::exp(b_val - a_val));
        let dldb = 1.0 / (1.0 + f64::exp(a_val - b_val));
        grad_w0 += *count as f64 * (dlda / w0 - dldb / (1.0 - w0));
        grad_c += *count as f64 * (dldb * (i_f64 / c - 1.0));
    }
    vec![grad_w0, grad_c]
}
