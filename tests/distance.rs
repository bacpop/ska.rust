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
        .args(["--threads", "2"])
        .assert()
        .stdout_eq_path(sandbox.file_string("merge.dist.stdout", TestDir::Correct));

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("distance")
        .arg(sandbox.file_string("merge.skf", TestDir::Input))
        .args(["-o", "dists.txt"])
        .assert()
        .success();
    sandbox.file_check("dists.txt", "merge.dist.stdout");

    // One difference
    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("distance")
        .arg(sandbox.file_string("merge_k41.skf", TestDir::Input))
        .arg("-v")
        .args(["--threads", "2"])
        .assert()
        .stdout_eq_path(sandbox.file_string("merge_k41.dist.stdout", TestDir::Correct));
}

// With filters
#[test]
fn dist_filter() {
    let sandbox = TestSetup::setup();

    // Has an S,G site, giving 0.5
    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("distance")
        .arg(sandbox.file_string("merge_k9.skf", TestDir::Input))
        .arg("--allow-ambiguous")
        .arg("-v")
        .args(["--threads", "2"])
        .assert()
        .stdout_eq_path(sandbox.file_string("merge_k9.dist.stdout", TestDir::Correct));

    // Ignoring the S,G gives 2.0
    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("distance")
        .arg(sandbox.file_string("merge_k9.skf", TestDir::Input))
        .arg("-v")
        .assert()
        .stdout_eq_path(sandbox.file_string("merge_k9_no_ambig.dist.stdout", TestDir::Correct));

    // Making only k-mers in both reduces mismatches to zero
    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("distance")
        .arg(sandbox.file_string("merge_k9.skf", TestDir::Input))
        .args(["--min-freq", "1"])
        .arg("-v")
        .assert()
        .stdout_eq_path(sandbox.file_string("merge_k9_min_freq.dist.stdout", TestDir::Correct));
}

#[test]
fn multisample_dists() {
    let sandbox = TestSetup::setup();

    // NB: Ensure dists of sample pair are the same as the above test
    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("build")
        .arg(sandbox.file_string("N_test_1.fa", TestDir::Input))
        .arg(sandbox.file_string("N_test_2.fa", TestDir::Input))
        .arg(sandbox.file_string("ambig_test_1.fa", TestDir::Input))
        .arg(sandbox.file_string("ambig_test_2.fa", TestDir::Input))
        .arg(sandbox.file_string("test_1.fa", TestDir::Input))
        .arg(sandbox.file_string("test_2.fa", TestDir::Input))
        .arg("-k")
        .arg("9")
        .arg("-o")
        .arg("multidist")
        .assert()
        .success();

    // Test with defaults (and that output order correct irrespective of threads)
    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("distance")
        .arg("multidist.skf")
        .arg("-v")
        .args(["--threads", "2"])
        .assert()
        .stdout_eq_path(sandbox.file_string("multidist.stdout", TestDir::Correct));

    // Test with min freq
    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("distance")
        .arg("multidist.skf")
        .arg("-v")
        .args(["--min-freq", "0.9"])
        .assert()
        .stdout_eq_path(sandbox.file_string("multidist.minfreq.stdout", TestDir::Correct));

    // Test with ambig bases
    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("distance")
        .arg("multidist.skf")
        .arg("-v")
        .arg("--allow-ambiguous")
        .assert()
        .stdout_eq_path(sandbox.file_string("multidist.ambig.stdout", TestDir::Correct));
}
