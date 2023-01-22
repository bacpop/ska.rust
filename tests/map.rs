use snapbox::cmd::{cargo_bin, Command};

pub mod common;
use crate::common::{TestDir, TestSetup};

// NB: to view output, uncomment the current_dir lines

#[test]
fn map_aln() {
    let sandbox = TestSetup::setup();

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("map")
        .arg(sandbox.file_string("test_ref.fa", TestDir::Input))
        .arg(sandbox.file_string("merge.skf", TestDir::Input))
        .args(&["-v", "--threads", "2"])
        .arg("-o")
        .arg("map.aln")
        .assert()
        .success();

    assert_eq!(true, sandbox.file_exists("map.aln"));

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("map")
        .arg(sandbox.file_string("test_ref.fa", TestDir::Input))
        .arg(sandbox.file_string("test_1.fa", TestDir::Input))
        .arg(sandbox.file_string("test_2.fa", TestDir::Input))
        .assert()
        .stdout_matches_path(sandbox.file_string("map_aln.stdout", TestDir::Correct));
}

#[test]
fn map_vcf() {
    let sandbox = TestSetup::setup();

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("map")
        .arg(sandbox.file_string("test_ref.fa", TestDir::Input))
        .arg(sandbox.file_string("merge.skf", TestDir::Input))
        .arg("-f")
        .arg("vcf")
        .arg("-o")
        .arg("map.vcf")
        .assert()
        .success();

    assert_eq!(true, sandbox.file_exists("map.vcf"));

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("map")
        .arg(sandbox.file_string("test_ref.fa", TestDir::Input))
        .arg(sandbox.file_string("test_1.fa", TestDir::Input))
        .arg(sandbox.file_string("test_2.fa", TestDir::Input))
        .arg("-f")
        .arg("vcf")
        .assert()
        .stdout_matches_path(sandbox.file_string("map_vcf.stdout", TestDir::Correct));
}
