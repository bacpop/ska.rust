use snapbox::cmd::{cargo_bin, Command};

mod common;
use crate::common::{TestDir, TestSetup};

#[test]
fn align_n() {
    let sandbox = TestSetup::setup();

    // Tests that N and n are skipped, so second variant doesn't appear in align
    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("build")
        .arg(sandbox.file_string("N_test_1.fa", TestDir::Input))
        .arg(sandbox.file_string("N_test_2.fa", TestDir::Input))
        .arg("-o")
        .arg("N_test")
        .assert()
        .success();

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("align")
        .arg("N_test.skf")
        .assert()
        .stdout_matches_path(sandbox.file_string("align_N.stdout", TestDir::Correct));
}

#[test]
fn map_n() {
    let sandbox = TestSetup::setup();

    // Tests that N and n are skipped, so second variant doesn't appear in map
    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("build")
        .arg(sandbox.file_string("N_test_1.fa", TestDir::Input))
        .arg(sandbox.file_string("N_test_2.fa", TestDir::Input))
        .arg("-o")
        .arg("N_test")
        .assert()
        .success();

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("map")
        .arg(sandbox.file_string("test_ref.fa", TestDir::Input))
        .arg("N_test.skf")
        .assert()
        .stdout_matches_path(sandbox.file_string("map_N.stdout", TestDir::Correct));
}

#[test]
fn repeats() {
    let sandbox = TestSetup::setup();

    // Tests that three repeats are correctly shown with IUPAC codes
    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("build")
        .arg("-k")
        .arg("9")
        .arg(sandbox.file_string("dup_test_1.fa", TestDir::Input))
        .arg(sandbox.file_string("dup_test_2.fa", TestDir::Input))
        .arg("-o")
        .arg("dup_ss")
        .arg("--single-strand")
        .assert()
        .success();

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("align")
        .arg("dup_ss.skf")
        .assert()
        .stdout_matches_path(sandbox.file_string("dup_ss.stdout", TestDir::Correct));

    // Also tests this is just a single variant
    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("weed")
        .arg("dup_ss.skf")
        .arg("--remove-const-sites")
        .args(&["--min-freq", "1"])
        .assert()
        .success();

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("nk")
        .arg("dup_ss.skf")
        .arg("--full-info")
        .assert()
        .stdout_matches_path(sandbox.file_string("dup_ss_nk.stdout", TestDir::Correct));

    // Tests that the RC of these IUPAC codes is correct
    // (the rc of the input k-mer is the canonical one)
    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("build")
        .arg("-k")
        .arg("9")
        .arg(sandbox.file_string("dup_test_1.fa", TestDir::Input))
        .arg(sandbox.file_string("dup_test_2.fa", TestDir::Input))
        .arg("-o")
        .arg("dup_rc")
        .assert()
        .success();

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("align")
        .arg("dup_rc.skf")
        .assert()
        .stdout_matches_path(sandbox.file_string("dup_rc.stdout", TestDir::Correct));
}
