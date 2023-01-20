use snapbox::cmd::{cargo_bin, Command};

mod common;
use crate::common::{TestDir, TestSetup};

// NB: to view output, uncomment the current_dir lines

#[test]
fn build_cli() {
    let sandbox = TestSetup::setup();
    // Create an rfile in the tmp dir
    let rfile_name = sandbox.create_rfile(false);

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("build")
        .arg("-f")
        .arg(rfile_name)
        .arg("-o")
        .arg("basic_build_opts")
        .args(&["-v", "--threads", "2", "-k", "31"])
        .assert()
        .success();

    assert_eq!(true, sandbox.file_exists("basic_build_opts.skf"));
}

#[test]
fn align_cli() {
    let sandbox = TestSetup::setup();

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("align")
        .arg(sandbox.file_string("test_1.fa", TestDir::Input))
        .arg(sandbox.file_string("test_2.fa", TestDir::Input))
        .arg("-o")
        .arg("basic.aln")
        .args(&["-v", "--threads", "2", "--const-sites", "--min-freq", "0"])
        .assert()
        .success();

    assert_eq!(true, sandbox.file_exists("basic.aln"));
}

#[test]
fn build_and_align() {
    let sandbox = TestSetup::setup();

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("build")
        .arg("-o")
        .arg("basic_build")
        .arg("-k")
        .arg("15")
        .arg(sandbox.file_string("test_1.fa", TestDir::Input))
        .arg(sandbox.file_string("test_2.fa", TestDir::Input))
        .assert()
        .success();

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("align")
        .arg("basic_build.skf")
        .arg("-o")
        .arg("basic.aln")
        .assert()
        .success();

    assert_eq!(
        true,
        sandbox.file_check("basic.aln", "basic1.aln")
            | sandbox.file_check("basic.aln", "basic2.aln")
    );
}

#[test]
fn basic_align() {
    let sandbox = TestSetup::setup();

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("align")
        .arg(sandbox.file_string("test_1.fa", TestDir::Input))
        .arg(sandbox.file_string("test_2.fa", TestDir::Input))
        .arg("-o")
        .arg("basic.aln")
        .assert()
        .success();

    assert_eq!(
        true,
        sandbox.file_check("basic.aln", "basic1.aln")
            | sandbox.file_check("basic.aln", "basic2.aln")
    );
}
