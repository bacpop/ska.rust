use snapbox::cmd::{cargo_bin, Command};

mod common;
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
