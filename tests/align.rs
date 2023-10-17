use snapbox::cmd::{cargo_bin, Command};

pub mod common;
use crate::common::*;

use hashbrown::HashSet;

#[cfg(test)]
use pretty_assertions::assert_eq;

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
        .args(&[
            "-v",
            "--threads",
            "2",
            "--filter",
            "no-filter",
            "--min-freq",
            "0",
        ])
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
        .arg("-v")
        .output()
        .unwrap()
        .stdout;

    let correct_aln = HashSet::from([vec!['A', 'T'], vec!['C', 'T']]);
    assert_eq!(var_hash(&fasta_align_out), correct_aln);
}

#[test]
fn long_kmers() {
    let sandbox = TestSetup::setup();

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("build")
        .arg("-o")
        .arg("build_k33")
        .arg("-k")
        .arg("33")
        .arg(sandbox.file_string("test_1.fa", TestDir::Input))
        .arg(sandbox.file_string("test_2.fa", TestDir::Input))
        .arg("-v")
        .assert()
        .success();

    let fasta_align_out_k33 = Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("align")
        .arg("build_k33.skf")
        .arg("-v")
        .output()
        .unwrap()
        .stdout;

    let correct_aln = HashSet::from([vec!['C', 'T'], vec!['T', 'A']]);
    assert_eq!(var_hash(&fasta_align_out_k33), correct_aln);

    // Check 128 bits used
    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("nk")
        .arg("build_k33.skf")
        .arg("-v")
        .assert()
        .stdout_matches_path(sandbox.file_string("k33.stdout", TestDir::Correct));

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("build")
        .arg("-o")
        .arg("build_k33")
        .arg("-k")
        .arg("65")
        .arg("-v")
        .arg(sandbox.file_string("test_1.fa", TestDir::Input))
        .arg(sandbox.file_string("test_2.fa", TestDir::Input))
        .assert()
        .failure();
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
        .arg(sandbox.file_string("merge_k9.skf", TestDir::Input))
        .arg("--filter")
        .arg("no-filter")
        .arg("--no-gap-only-sites") // adding with no filter produces a warning
        .arg("-v")
        .output()
        .unwrap()
        .stdout;

    let full_length = 38; // could also get 38 from ska nk
    let lengths = aln_length(&unfilt_align_out);
    for length in lengths {
        assert_eq!(full_length, length);
    }

    let noambig_align_out = Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("align")
        .arg(sandbox.file_string("merge_k9.skf", TestDir::Input))
        .arg("--filter")
        .arg("no-ambig")
        .arg("-v")
        .output()
        .unwrap()
        .stdout;

    let lengths = aln_length(&noambig_align_out);
    for length in lengths {
        assert_eq!(full_length - 1, length);
    }

    let const_filt_align_out = Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("align")
        .arg(sandbox.file_string("merge_k9.skf", TestDir::Input))
        .arg("--filter")
        .arg("no-const")
        .arg("-v")
        .output()
        .unwrap()
        .stdout;

    let correct_aln = HashSet::from([vec!['T', 'A'], vec!['C', 'T'], vec!['S', 'G']]);
    assert_eq!(var_hash(&const_filt_align_out), correct_aln);

    let no_ambig_or_const_align_out = Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("align")
        .arg(sandbox.file_string("merge_k9.skf", TestDir::Input))
        .arg("--filter")
        .arg("no-ambig-or-const")
        .arg("-v")
        .output()
        .unwrap()
        .stdout;

    let correct_aln = HashSet::from([vec!['T', 'A'], vec!['C', 'T']]);
    assert_eq!(var_hash(&no_ambig_or_const_align_out), correct_aln);

    let ambig_filter_align_out = Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("align")
        .arg(sandbox.file_string("merge_k9.skf", TestDir::Input))
        .arg("--filter")
        .arg("no-const")
        .arg("--ambig-mask")
        .arg("-v")
        .output()
        .unwrap()
        .stdout;

    let correct_aln = HashSet::from([vec!['T', 'A'], vec!['C', 'T'], vec!['N', 'G']]);
    assert_eq!(var_hash(&ambig_filter_align_out), correct_aln);

    // Output everything by using min-freq 0, check for behaviour on gap only sites
    let unfilt_align_out = Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("align")
        .arg(sandbox.file_string("merge_k9.skf", TestDir::Input))
        .arg("--filter")
        .arg("no-const")
        .arg("--min-freq")
        .arg("0")
        .arg("-v")
        .output()
        .unwrap()
        .stdout;

    let full_length = 33;
    let lengths = aln_length(&unfilt_align_out);
    for length in lengths {
        assert_eq!(full_length, length);
    }

    let filt_align_out = Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("align")
        .arg(sandbox.file_string("merge_k9.skf", TestDir::Input))
        .arg("--filter")
        .arg("no-const")
        .arg("--min-freq")
        .arg("0")
        .arg("--no-gap-only-sites")
        .arg("-v")
        .output()
        .unwrap()
        .stdout;

    let full_length = 3;
    let lengths = aln_length(&filt_align_out);
    for length in lengths {
        assert_eq!(full_length, length);
    }

    let unfilt_align_out = Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("align")
        .arg(sandbox.file_string("merge_k9.skf", TestDir::Input))
        .arg("--filter")
        .arg("no-ambig-or-const")
        .arg("--min-freq")
        .arg("0")
        .arg("-v")
        .output()
        .unwrap()
        .stdout;

    let full_length = 32;
    let lengths = aln_length(&unfilt_align_out);
    for length in lengths {
        assert_eq!(full_length, length);
    }

    let filt_align_out = Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("align")
        .arg(sandbox.file_string("merge_k9.skf", TestDir::Input))
        .arg("--filter")
        .arg("no-ambig-or-const")
        .arg("--min-freq")
        .arg("0")
        .arg("-v")
        .arg("--no-gap-only-sites")
        .output()
        .unwrap()
        .stdout;

    let full_length = 2;
    let lengths = aln_length(&filt_align_out);
    for length in lengths {
        assert_eq!(full_length, length);
    }
}

#[test]
fn parallel_align() {
    let sandbox = TestSetup::setup();
    // Needs at least ten samples per threads, so make lots of copies
    // (confirmed with -v that this actually uses the multi-thread algorithm)
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
        .args(&["-v", "--threads", "4", "-k", "15"])
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
