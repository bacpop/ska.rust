use snapbox::cmd::{cargo_bin, Command};

pub mod common;
use crate::common::*;

use hashbrown::HashSet;

// NB: to view output, uncomment the current_dir lines

#[test]
fn build_cli() {
    let sandbox = TestSetup::setup();
    // Create an rfile in the tmp dir
    let rfile_name = sandbox.create_rfile(&"test", FxType::Fasta);

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
        .args(&["-v", "--threads", "2", "--filter", "no-filter", "--min-freq", "0"])
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

    let fasta_align_out = Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("align")
        .arg("basic_build.skf")
        .output()
        .unwrap()
        .stdout;

    let correct_aln = HashSet::from([vec!['A', 'T'], vec!['C', 'T']]);
    assert_eq!(var_hash(&fasta_align_out), correct_aln);
}

#[test]
fn basic_align() {
    let sandbox = TestSetup::setup();

    let fasta_align_out = Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("align")
        .arg(sandbox.file_string("test_1.fa", TestDir::Input))
        .arg(sandbox.file_string("test_2.fa", TestDir::Input))
        .output()
        .unwrap()
        .stdout;

    let correct_aln = HashSet::from([vec!['A', 'T'], vec!['C', 'T']]);
    assert_eq!(var_hash(&fasta_align_out), correct_aln);
}

#[test]
fn filters() {
    let sandbox = TestSetup::setup();

    // With k=9 there is a repeated k-mer
    let unfilt_align_out = Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("align")
        .arg(sandbox.file_string("merge_k9.fa", TestDir::Input))
        .arg("--filter")
        .arg("no-filter")
        .output()
        .unwrap()
        .stdout;

    let lengths = aln_length(&unfilt_align_out);
    for length in lengths {
        assert_eq!(39, length); // could also get 39 from ska nk
    }

    let const_filt_align_out = Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("align")
        .arg(sandbox.file_string("merge_k9.skf", TestDir::Input))
        .arg("--filter")
        .arg("no-const")
        .output()
        .unwrap()
        .stdout;


    let mut correct_aln = HashSet::from([vec!['T', 'A'], vec!['C', 'T'], vec!['S', 'G']]);
    assert_eq!(var_hash(&const_filt_align_out), correct_aln);

    let no_ambig_align_out = Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("align")
        .arg(sandbox.file_string("merge_k9.skf", TestDir::Input))
        .arg("--filter")
        .arg("no-ambig-or-const")
        .output()
        .unwrap()
        .stdout;

    correct_aln.remove(&vec!['S', 'G']);
    assert_eq!(var_hash(&no_ambig_align_out), correct_aln);
}

#[test]
fn parallel_align() {
    let sandbox = TestSetup::setup();
    // Needs at least ten samples per threads, so make lots of copies
    // (confimed with -v that this actually uses the multi-thread algorithm)
    let rfile_name = sandbox.create_par_rfile();

    // Serial alignment
    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("build")
        .arg("-f")
        .arg(rfile_name)
        .arg("-o")
        .arg("serial_build")
        .args(&["-v", "--threads", "1", "-k", "15"])
        .assert()
        .success();

    let serial_align_out = Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("align")
        .arg("serial_build.skf")
        .output()
        .unwrap()
        .stdout;

    // Parallel alignment algorithm used
    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("build")
        .arg("-f")
        .arg(rfile_name)
        .arg("-o")
        .arg("parallel_build")
        .args(&["-v", "--threads", "2", "-k", "15"])
        .assert()
        .success();

    let parallel_align_out = Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("align")
        .arg("parallel_build.skf")
        .output()
        .unwrap()
        .stdout;

    assert_eq!(var_hash(&serial_align_out), var_hash(&parallel_align_out));
}
