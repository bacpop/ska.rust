use snapbox::cmd::{cargo_bin, Command};

#[cfg(test)]
use pretty_assertions::{assert_eq};

pub mod common;
use crate::common::{TestDir, TestSetup};

// NB: to view output, uncomment the current_dir lines

#[test]
fn merge_delete() {
    let sandbox = TestSetup::setup();

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("build")
        .arg(sandbox.file_string("test_1.fa", TestDir::Input))
        .arg("-o")
        .arg("test_1")
        .assert()
        .success();

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("build")
        .arg(sandbox.file_string("test_2.fa", TestDir::Input))
        .arg("-o")
        .arg("test_2")
        .assert()
        .success();

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("merge")
        .arg("test_1.skf")
        .arg("test_2.skf")
        .arg("-o")
        .arg("merge")
        .assert()
        .success();

    assert_eq!(true, sandbox.file_exists("merge.skf"));

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("nk")
        .arg("merge.skf")
        .assert()
        .stdout_matches_path(sandbox.file_string("merge_nk.stdout", TestDir::Correct));

    // Try removing non-existent sample
    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("delete")
        .arg("test_3")
        .arg("merge.skf")
        .assert()
        .failure();

    // delete a sample then check nk same as original
    let test1_nk = Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("nk")
        .arg("test_1.skf")
        .output()
        .unwrap();

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("delete")
        .arg("merge.skf")
        .arg("test_2")
        .assert()
        .success();

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("nk")
        .arg("merge.skf")
        .assert()
        .stdout_eq(test1_nk.stdout);
}

#[test]
fn merge_delete_u128() {
    let sandbox = TestSetup::setup();

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("build")
        .arg(sandbox.file_string("test_1.fa", TestDir::Input))
        .arg("-o")
        .arg("test_1")
        .arg("-k")
        .arg("41")
        .assert()
        .success();

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("build")
        .arg(sandbox.file_string("test_2.fa", TestDir::Input))
        .arg("-o")
        .arg("test_2")
        .arg("-k")
        .arg("41")
        .assert()
        .success();

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("merge")
        .arg("test_1.skf")
        .arg("test_2.skf")
        .arg("-o")
        .arg("merge")
        .arg("-v")
        .assert()
        .success();

    assert_eq!(true, sandbox.file_exists("merge.skf"));

    // Try removing non-existent sample
    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("delete")
        .arg("test_3")
        .arg("merge.skf")
        .arg("-v")
        .assert()
        .failure();

    // delete a sample then check nk same as original
    let test1_nk = Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("nk")
        .arg("test_1.skf")
        .output()
        .unwrap();

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("delete")
        .arg("merge.skf")
        .arg("test_2")
        .arg("-v")
        .assert()
        .success();

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("nk")
        .arg("merge.skf")
        .assert()
        .stdout_eq(test1_nk.stdout);
}

#[test]
fn weed() {
    let sandbox = TestSetup::setup();

    Command::new("cp")
        .current_dir(sandbox.get_wd())
        .arg(sandbox.file_string("merge.skf", TestDir::Input))
        .arg("merge.skf")
        .assert()
        .success();

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("weed")
        .arg("merge.skf")
        .arg(sandbox.file_string("weed.fa", TestDir::Input))
        .arg("-v")
        .assert()
        .success();

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("align")
        .arg("merge.skf")
        .assert()
        .stdout_eq_path(sandbox.file_string("weed_align.stdout", TestDir::Correct));

    // With const sites/filter and nk full info
    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("weed")
        .arg("merge.skf")
        .arg("--filter")
        .arg("no-const")
        .args(&["--min-freq", "1"])
        .assert()
        .success();

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("nk")
        .arg("merge.skf")
        .arg("--full-info")
        .assert()
        .stdout_matches_path(sandbox.file_string("weed_nk.stdout", TestDir::Correct));

    // With longer k-mers
    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("build")
        .arg("-o")
        .arg("build_k41")
        .arg("-k")
        .arg("41")
        .arg(sandbox.file_string("test_1.fa", TestDir::Input))
        .arg(sandbox.file_string("test_2.fa", TestDir::Input))
        .arg("-v")
        .assert()
        .success();

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("weed")
        .arg("build_k41.skf")
        .arg("--filter")
        .arg("no-ambig-or-const")
        .arg("-v")
        .assert()
        .success();

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("nk")
        .arg("build_k41.skf")
        .arg("--full-info")
        .arg("-v")
        .assert()
        .stdout_matches_path(sandbox.file_string("weed_nk_k41.stdout", TestDir::Correct));
}
