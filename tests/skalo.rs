use snapbox::cmd::{cargo_bin, Command};

#[cfg(test)]
pub mod common;
use crate::common::*;

// NB: to view output, uncomment the current_dir lines

// reference-free SNP/indel calling with positioning on a reference genome
#[test]
fn ska_lo() {
    let sandbox = TestSetup::setup();

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("lo")
        .arg("-r")
        .arg(sandbox.file_string("test_skalo_reference.fas", TestDir::Input))
        .arg(sandbox.file_string("test_skalo.skf", TestDir::Input))
        .arg("test_skalo")
        .assert()
        .success();

    sandbox.file_check("test_skalo_snps.fas", "test_skalo_snps.fas");

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("lo")
        .arg(sandbox.file_string("test_skalo_indel.skf", TestDir::Input))
        .arg("test_skalo")
        .assert()
        .success();

    sandbox.file_check("test_skalo_indels.vcf", "test_skalo_indels.vcf");
}
