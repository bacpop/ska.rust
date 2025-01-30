use snapbox::cmd::{cargo_bin, Command};

#[cfg(test)]
use pretty_assertions::assert_eq;

pub mod common;
use crate::common::*;

// NB: to view output, uncomment the current_dir lines

#[test]
fn ska_lo() {
    let sandbox = TestSetup::setup();

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("lo")
        .arg("-i")
        .arg(sandbox.file_string("test_skalo.skf", TestDir::Input))
        .arg("-o")
        .arg("test_skalo")
        .assert()
        .success();

    assert_eq!(true, sandbox.file_exists("test_skalo_snps.fas"));
}
