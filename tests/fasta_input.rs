use snapbox::cmd::{cargo_bin, Command};

use hashbrown::HashSet;

pub mod common;
use crate::common::{var_hash, TestDir, TestSetup};

#[cfg(test)]
use pretty_assertions::assert_eq;

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
        .arg("N_test.skf") // test no .skf.skf extension
        .assert()
        .success();

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("align")
        .arg("N_test.skf")
        .assert()
        .stdout_eq_path(sandbox.file_string("align_N.stdout", TestDir::Correct));
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
        .arg("-k")
        .arg("11")
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
        .stdout_eq_path(sandbox.file_string("map_N.stdout", TestDir::Correct));
}

#[test]
fn rev_comp() {
    let sandbox = TestSetup::setup();

    // Test an RC sequence gives same alignment out
    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("build")
        .arg("-o")
        .arg("fwd_build")
        .arg("-k")
        .arg("15")
        .arg(sandbox.file_string("test_1.fa", TestDir::Input))
        .arg(sandbox.file_string("test_2.fa", TestDir::Input))
        .assert()
        .success();

    let no_rc_aln = Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("align")
        .arg("fwd_build.skf")
        .output()
        .unwrap()
        .stdout;

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("build")
        .arg("-o")
        .arg("fwd_build")
        .arg("-k")
        .arg("15")
        .arg(sandbox.file_string("test_1.fa", TestDir::Input))
        .arg(sandbox.file_string("test_2_rc.fa", TestDir::Input))
        .assert()
        .success();

    let rc_aln = Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("align")
        .arg("fwd_build.skf")
        .output()
        .unwrap()
        .stdout;

    assert_eq!(var_hash(&no_rc_aln), var_hash(&rc_aln));

    // Now with rc off
    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("build")
        .arg("-o")
        .arg("fwd_build")
        .arg("-k")
        .arg("15")
        .arg(sandbox.file_string("test_1.fa", TestDir::Input))
        .arg(sandbox.file_string("test_2_rc.fa", TestDir::Input))
        .arg("--single-strand")
        .assert()
        .success();

    let ss_aln = Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("align")
        .arg("fwd_build.skf")
        .output()
        .unwrap()
        .stdout;

    assert_eq!(var_hash(&ss_aln).is_empty(), true);

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("build")
        .arg("-o")
        .arg("build_k33")
        .arg("-k")
        .arg("33")
        .arg("--single-strand")
        .arg(sandbox.file_string("test_1.fa", TestDir::Input))
        .arg(sandbox.file_string("test_2.fa", TestDir::Input))
        .assert()
        .success();

    let fasta_align_out_k33 = Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("align")
        .arg("build_k33.skf")
        .output()
        .unwrap()
        .stdout;

    let correct_aln = HashSet::from([vec!['T', 'A'], vec!['G', 'A']]);
    assert_eq!(var_hash(&fasta_align_out_k33), correct_aln);
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
        .stdout_eq_path(sandbox.file_string("dup_ss.stdout", TestDir::Correct));

    // Also tests this is just a single variant
    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("weed")
        .arg("dup_ss.skf")
        .arg("--filter")
        .arg("no-const")
        .args(["--min-freq", "1"])
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
        .stdout_eq_path(sandbox.file_string("dup_rc.stdout", TestDir::Correct));
}

#[test]
fn palindromes() {
    let sandbox = TestSetup::setup();

    // Check that palindrome k-mers (self-rc) are correctly ambiguous
    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("build")
        .arg("-k")
        .arg("15")
        .arg(sandbox.file_string("palindrome_1.fa", TestDir::Input))
        .arg(sandbox.file_string("palindrome_2.fa", TestDir::Input))
        .arg("-o")
        .arg("otto")
        .arg("-v")
        .assert()
        .success();

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("align")
        .arg("otto.skf")
        .args(["--filter", "no-filter"])
        .assert()
        .stdout_eq_path(sandbox.file_string("palindrome.stdout", TestDir::Correct));

    // Single strand forces orientation so no ambiguity
    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("build")
        .arg("-k")
        .arg("15")
        .arg(sandbox.file_string("palindrome_1.fa", TestDir::Input))
        .arg(sandbox.file_string("palindrome_2.fa", TestDir::Input))
        .arg("-o")
        .arg("otan")
        .arg("--single-strand")
        .assert()
        .success();

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("align")
        .arg("otan.skf")
        .assert()
        .stdout_eq_path(sandbox.file_string("palindrome_norc.stdout", TestDir::Correct));

    // Check that palindrome k-mers (self-rc) are correctly ambiguous, and
    // work correctly when multiple copies are present
    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("build")
        .arg("-k")
        .arg("15")
        .arg(sandbox.file_string("palindrome_reps_1.fa", TestDir::Input))
        .arg(sandbox.file_string("palindrome_reps_2.fa", TestDir::Input))
        .arg("-o")
        .arg("ottootto")
        .arg("-v")
        .assert()
        .success();

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("align")
        .arg("ottootto.skf")
        .args(["--filter", "no-filter"])
        .assert()
        .stdout_eq_path(sandbox.file_string("palindrome_reps.stdout", TestDir::Correct));
}
