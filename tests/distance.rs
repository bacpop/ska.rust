use snapbox::cmd::{cargo_bin, Command};

pub mod common;
use crate::common::*;

// NB: to view output, uncomment the current_dir lines

// With two samples
#[test]
fn basic_dists() {
    let sandbox = TestSetup::setup();

    // Two differences, so dist 2.0
    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("distance")
        .arg(sandbox.file_string("merge.skf", TestDir::Input))
        .arg("-v")
        .args(&["--threads", "2"])
        .assert()
        .stdout_eq_path(sandbox.file_string("merge.dist.stdout", TestDir::Correct));

    // One difference
    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("distance")
        .arg(sandbox.file_string("merge_k41.skf", TestDir::Input))
        .arg("-v")
        .args(&["--threads", "2"])
        .assert()
        .stdout_eq_path(sandbox.file_string("merge_k41.dist.stdout", TestDir::Correct));
}

#[test]
fn dist_filter() {
    let sandbox = TestSetup::setup();

    // Has an S,G site, giving 0.5
    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("distance")
        .arg(sandbox.file_string("merge_k9.skf", TestDir::Input))
        .arg("-v")
        .args(&["--threads", "2"])
        .assert()
        .stdout_eq_path(sandbox.file_string("merge_k9.dist.stdout", TestDir::Correct));

    // Ignoring the S,G gives 2.0
    // TODO mismatches seem wrong here
    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("distance")
        .arg(sandbox.file_string("merge_k9.skf", TestDir::Input))
        .arg("--filter-ambiguous")
        .arg("-v")
        .assert()
        .stdout_eq_path(sandbox.file_string("merge_k9_no_ambig.dist.stdout", TestDir::Correct));

    // Making only k-mers in both reduces mismatches to zero
    // TODO why does this change to 2.0?
    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("distance")
        .arg(sandbox.file_string("merge_k9.skf", TestDir::Input))
        .args(&["--min-freq", "1"])
        .arg("-v")
        .assert()
        .stdout_eq_path(sandbox.file_string("merge_k9_min_freq.dist.stdout", TestDir::Correct));
}

// TODO need a multisample test