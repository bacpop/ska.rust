use std::path::Path;
use predicates::prelude::*;

use snapbox::cmd::{Command, cargo_bin};

// Creates correct path for input/output files
static FILE_IN: &'static str = "tests/files_in";
static FILE_OUT: &'static str = "tests/files_out";
static FILE_TEST: &'static str = "tests/files_test";

#[macro_export]
macro_rules! file_in {
    ($e:expr) => {&format!("{}/{}", FILE_IN, $e)}
}

#[macro_export]
macro_rules! file_out {
    ($e:expr) => {&format!("{}/{}", FILE_OUT, $e)}
}

#[macro_export]
macro_rules! file_test {
    ($e:expr) => {&format!("{}/{}", FILE_TEST, $e)}
}

#[test]
fn build_and_align() {
    Command::new(cargo_bin("ska"))
        .arg("build")
        .arg("-o")
        .arg(file_out!("basic_build"))
        .arg("-k")
        .arg("17")
        .arg(file_in!("test_1.fa"))
        .arg(file_in!("test_2.fa"))
        .assert()
        .success();

    Command::new(cargo_bin("ska"))
        .arg("align")
        .arg(file_out!("basic_build.skf"))
        .arg("-o")
        .arg(file_out!("basic.aln"))
        .assert()
        .success();

    let predicate_file = predicate::path::eq_file(Path::new(file_out!("basic.aln")));
    assert_eq!(true, predicate_file.eval(Path::new(file_test!("basic1.aln"))) |
                     predicate_file.eval(Path::new(file_test!("basic2.aln"))));
}

#[test]
fn basic_align() {
    Command::new(cargo_bin("ska"))
        .arg("align")
        .arg(file_in!("test_1.fa"))
        .arg(file_in!("test_2.fa"))
        .arg("-o")
        .arg(file_out!("basic.aln"))
        .assert()
        .success();

    let predicate_file = predicate::path::eq_file(Path::new(file_out!("basic.aln")));
    assert_eq!(true, predicate_file.eval(Path::new(file_test!("basic1.aln"))) |
                     predicate_file.eval(Path::new(file_test!("basic2.aln"))));
}

