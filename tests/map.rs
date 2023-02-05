use snapbox::cmd::{cargo_bin, Command};

pub mod common;
use crate::common::*;

// NB: to view output, uncomment the current_dir lines

#[test]
fn map_aln() {
    let sandbox = TestSetup::setup();

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("map")
        .arg(sandbox.file_string("test_ref.fa", TestDir::Input))
        .arg(sandbox.file_string("test_1.fa", TestDir::Input))
        .arg(sandbox.file_string("test_2.fa", TestDir::Input))
        .arg("-o")
        .arg("map.aln")
        .assert()
        .success();

    assert_eq!(true, sandbox.file_exists("map.aln"));

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("map")
        .arg(sandbox.file_string("test_ref.fa", TestDir::Input))
        .arg(sandbox.file_string("merge.skf", TestDir::Input))
        .args(&["-v", "--threads", "2"])
        .assert()
        .stdout_matches_path(sandbox.file_string("map_aln.stdout", TestDir::Correct));

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("map")
        .arg(sandbox.file_string("test_ref.fa", TestDir::Input))
        .arg(sandbox.file_string("merge_k9.skf", TestDir::Input))
        .args(&["-v", "--threads", "2"])
        .assert()
        .stdout_matches_path(sandbox.file_string("map_aln_k9.stdout", TestDir::Correct));
}

#[test]
fn map_vcf() {
    let sandbox = TestSetup::setup();

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("map")
        .arg(sandbox.file_string("test_ref.fa", TestDir::Input))
        .arg(sandbox.file_string("test_1.fa", TestDir::Input))
        .arg(sandbox.file_string("test_2.fa", TestDir::Input))
        .arg("-o")
        .arg("map.vcf")
        .arg("-f")
        .arg("vcf")
        .assert()
        .success();

    assert_eq!(true, sandbox.file_exists("map.vcf"));

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("map")
        .arg(sandbox.file_string("test_ref.fa", TestDir::Input))
        .arg(sandbox.file_string("merge.skf", TestDir::Input))
        .arg("-f")
        .arg("vcf")
        .assert()
        .stdout_matches_path(sandbox.file_string("map_vcf.stdout", TestDir::Correct));
}

#[test]
fn map_rev_comp() {
    let sandbox = TestSetup::setup();

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("build")
        .arg("-o")
        .arg("rc_build")
        .arg("-k")
        .arg("9")
        .arg(sandbox.file_string("test_1.fa", TestDir::Input))
        .arg(sandbox.file_string("test_2_rc.fa", TestDir::Input))
        .assert()
        .success();

    let rc_map = Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("map")
        .arg(sandbox.file_string("test_ref.fa", TestDir::Input))
        .arg("fwd_build.skf")
        .output()
        .unwrap()
        .stdout;

    let fwd_map = Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("map")
        .arg(sandbox.file_string("test_ref.fa", TestDir::Input))
        .arg(sandbox.file_string("merge.skf", TestDir::Input))
        .output()
        .unwrap()
        .stdout;

    cmp_map_aln(&rc_map, &fwd_map);

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("build")
        .arg("-o")
        .arg("ss_map")
        .arg("-k")
        .arg("9")
        .arg("--single-strand")
        .arg(sandbox.file_string("test_1.fa", TestDir::Input))
        .arg(sandbox.file_string("test_2_rc.fa", TestDir::Input))
        .assert()
        .success();
}
