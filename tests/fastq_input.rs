use snapbox::cmd::{cargo_bin, Command};

use hashbrown::HashSet;

#[cfg(test)]
use pretty_assertions::assert_eq;

pub mod common;
use crate::common::*;

// Uses perfect reads with the same fasta input as merge.skf
#[test]
fn align_fastq() {
    let sandbox = TestSetup::setup();
    let rfile_name = sandbox.create_rfile("test", FxType::Fastq);

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("build")
        .arg("-f")
        .arg(rfile_name)
        .arg("-o")
        .arg("reads")
        .args(["--min-count", "2", "-v", "-k", "9", "--min-qual", "2"])
        .assert()
        .success();

    let fastq_align_out = Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("align")
        .arg("reads.skf")
        .output()
        .unwrap()
        .stdout;

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("build")
        .args(["-k", "9"])
        .arg(sandbox.file_string("test_1.fa", TestDir::Input))
        .arg(sandbox.file_string("test_2.fa", TestDir::Input))
        .args(["-o", "fasta_k9"])
        .output()
        .unwrap();

    let fasta_align_out = Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("align")
        .arg("fasta_k9.skf")
        .output()
        .unwrap()
        .stdout;

    assert_eq!(var_hash(&fastq_align_out), var_hash(&fasta_align_out));
}

// Case with errors in middle base, countmin should use both split k-mer
// and base to filter
#[test]
fn count_check() {
    let sandbox = TestSetup::setup();
    let rfile_name = sandbox.create_rfile("test_count", FxType::Fastq);

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("build")
        .arg("-f")
        .arg(rfile_name)
        .arg("-o")
        .arg("reads_k7_c1")
        .args(["--min-count", "1", "-v", "-k", "7"])
        .assert()
        .success();

    let fastq_align_c1_out = Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("align")
        .arg("reads_k7_c1.skf")
        .output()
        .unwrap()
        .stdout;

    let correct_aln = HashSet::from([vec!['C', 'W']]);
    assert_eq!(var_hash(&fastq_align_c1_out), correct_aln);

    // In sample two there are three split k-mers supporting T, one supporting A
    // So above see a W, here the A gets filtered and we just see a T
    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("build")
        .arg("-f")
        .arg(rfile_name)
        .arg("-o")
        .arg("reads_k7_c3")
        .args(["--min-count", "3", "-v", "-k", "7"])
        .assert()
        .success();

    let fastq_align_c3_out = Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("align")
        .arg("reads_k7_c3.skf")
        .output()
        .unwrap()
        .stdout;

    let correct_aln = HashSet::from([vec!['C', 'T']]);
    assert_eq!(var_hash(&fastq_align_c3_out), correct_aln);
}

// Case with errors in middle base, countmin should use both split k-mer
// and base to filter
#[test]
fn count_check_long() {
    let sandbox = TestSetup::setup();
    let rfile_name = sandbox.create_rfile("test_long", FxType::Fastq);

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("build")
        .arg("-f")
        .arg(rfile_name)
        .arg("-o")
        .arg("reads_k63_c1")
        .args(["--min-count", "1", "-v", "-k", "63"])
        .assert()
        .success();

    let fastq_align_c1_out = Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("align")
        .arg("reads_k63_c1.skf")
        .arg("-v")
        .output()
        .unwrap()
        .stdout;

    let correct_aln = HashSet::from([vec!['G', 'M']]);
    assert_eq!(var_hash(&fastq_align_c1_out), correct_aln);

    // In sample two there are three split k-mers supporting T, one supporting A
    // So above see a W, here the A gets filtered and we just see a T
    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("build")
        .arg("-f")
        .arg(rfile_name)
        .arg("-o")
        .arg("reads_k63_c3")
        .args(["--min-count", "3", "-v", "-k", "63"])
        .assert()
        .success();

    let fastq_align_c3_out = Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("align")
        .arg("reads_k63_c3.skf")
        .arg("-v")
        .output()
        .unwrap()
        .stdout;

    let correct_aln = HashSet::from([vec!['G', 'A']]);
    assert_eq!(var_hash(&fastq_align_c3_out), correct_aln);

    // Now ignoring the rc
    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("build")
        .arg("-f")
        .arg(rfile_name)
        .arg("-o")
        .arg("reads_k63_c3_ss")
        .arg("--single-strand")
        .args(["--min-count", "2", "-v", "-k", "63"])
        .assert()
        .success();

    let fastq_align_c3_out_ss = Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("align")
        .arg("reads_k63_c3_ss.skf")
        .arg("-v")
        .output()
        .unwrap()
        .stdout;

    assert_eq!(
        var_hash(&fastq_align_c3_out),
        var_hash(&fastq_align_c3_out_ss)
    );
}

// Uses perfect reads with the same fasta input as merge.skf
#[test]
fn map_fastq() {
    let sandbox = TestSetup::setup();
    let rfile_name = sandbox.create_rfile("test", FxType::Fastq);

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("build")
        .arg("-f")
        .arg(rfile_name)
        .arg("-o")
        .arg("reads")
        .args(["--min-count", "1", "-v", "-k", "9", "--min-qual", "2"])
        .assert()
        .success();

    let fastq_map_out = Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("map")
        .arg(sandbox.file_string("test_ref.fa", TestDir::Input))
        .arg("reads.skf")
        .output()
        .unwrap()
        .stdout;

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("build")
        .arg(sandbox.file_string("test_1.fa", TestDir::Input))
        .arg(sandbox.file_string("test_2.fa", TestDir::Input))
        .arg("-o")
        .arg("assemblies")
        .args(["-v", "-k", "9"])
        .assert()
        .success();

    let fasta_map_out = Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("map")
        .arg(sandbox.file_string("test_ref.fa", TestDir::Input))
        .arg("assemblies.skf")
        .output()
        .unwrap()
        .stdout;

    assert_eq!(
        String::from_utf8(fastq_map_out).unwrap(),
        String::from_utf8(fasta_map_out).unwrap()
    );

    let fastq_map_out_vcf = Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("map")
        .arg(sandbox.file_string("test_ref.fa", TestDir::Input))
        .arg("reads.skf")
        .args(["-f", "vcf"])
        .output()
        .unwrap()
        .stdout;

    let fasta_map_out_vcf = Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("map")
        .arg(sandbox.file_string("test_ref.fa", TestDir::Input))
        .arg("assemblies.skf")
        .args(["-f", "vcf"])
        .output()
        .unwrap()
        .stdout;

    assert_eq!(
        String::from_utf8(fastq_map_out_vcf),
        String::from_utf8(fasta_map_out_vcf)
    );
}

// Add errors and low quality scores to split k-mers around
// CACT    TTAA    C,T
// zgrep 'TTAA.AGTG' test_1_fwd.fastq.gz
// zgrep 'CACT.TTAA' test_1_rev.fastq.gz
// test_error has base errors
// test_quality has many qual '0' in these k-mers
// test_quality_base has qual '0' at the middle base and '!' at one flanking base
#[test]
fn error_fastq() {
    let sandbox = TestSetup::setup();

    // Without errors
    let rfile_name = sandbox.create_rfile("test", FxType::Fastq);
    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("build")
        .arg("-f")
        .arg(rfile_name)
        .arg("-o")
        .arg("reads")
        .args(["--min-count", "3", "-v", "-k", "9", "--min-qual", "2"])
        .assert()
        .success();

    let fastq_align_out_all = Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("align")
        .arg("reads.skf")
        .output()
        .unwrap()
        .stdout;
    let mut all_hash = var_hash(&fastq_align_out_all);

    // With no quality filtering
    let rfile_name = sandbox.create_rfile("test_quality", FxType::Fastq);
    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("build")
        .arg("-f")
        .arg(rfile_name)
        .arg("-o")
        .arg("reads")
        .args([
            "--min-count",
            "5",
            "-v",
            "-k",
            "9",
            "--qual-filter",
            "no-filter",
        ])
        .assert()
        .success();

    let fastq_align_out_quality_nofilter = Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("align")
        .arg("reads.skf")
        .output()
        .unwrap()
        .stdout;

    assert_eq!(var_hash(&fastq_align_out_quality_nofilter), all_hash);

    // With only quality filtering the middle base
    let rfile_name = sandbox.create_rfile("test_quality_base", FxType::Fastq);
    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("build")
        .arg("-f")
        .arg(rfile_name)
        .arg("-o")
        .arg("reads")
        .args([
            "--min-count",
            "5",
            "-v",
            "-k",
            "9",
            "--qual-filter",
            "middle",
            "--min-qual",
            "5",
        ])
        .assert()
        .success();

    let fastq_align_out_quality_base = Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("align")
        .arg("reads.skf")
        .output()
        .unwrap()
        .stdout;

    assert_eq!(var_hash(&fastq_align_out_quality_base), all_hash);

    // With errors
    all_hash.remove(&vec!['C', 'T']);
    let rfile_name = sandbox.create_rfile("test_error", FxType::Fastq);
    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("build")
        .arg("-f")
        .arg(rfile_name)
        .arg("-o")
        .arg("reads")
        .args(["--min-count", "5", "-v", "-k", "9", "--min-qual", "2"])
        .assert()
        .success();

    let fastq_align_out_error = Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("align")
        .arg("reads.skf")
        .output()
        .unwrap()
        .stdout;

    assert_eq!(var_hash(&fastq_align_out_error), all_hash);

    // With low quality score
    let rfile_name = sandbox.create_rfile("test_quality", FxType::Fastq);
    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("build")
        .arg("-f")
        .arg(rfile_name)
        .arg("-o")
        .arg("reads")
        .args(["--min-count", "5", "-v", "-k", "9", "--min-qual", "30"])
        .assert()
        .success();

    let fastq_align_out_quality = Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("align")
        .arg("reads.skf")
        .output()
        .unwrap()
        .stdout;

    assert_eq!(var_hash(&fastq_align_out_quality), all_hash);

    // With low quality score in flanking region
    let rfile_name = sandbox.create_rfile("test_quality_base", FxType::Fastq);
    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("build")
        .arg("-f")
        .arg(rfile_name)
        .arg("-o")
        .arg("reads")
        .args([
            "--min-count",
            "5",
            "-v",
            "-k",
            "9",
            "--min-qual",
            "5",
            "--qual-filter",
            "strict",
        ])
        .assert()
        .success();

    let fastq_align_out_quality = Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("align")
        .arg("reads.skf")
        .output()
        .unwrap()
        .stdout;

    assert_eq!(var_hash(&fastq_align_out_quality), all_hash);

    // With low quality score in middle base region
    let rfile_name = sandbox.create_rfile("test_quality_base", FxType::Fastq);
    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("build")
        .arg("-f")
        .arg(rfile_name)
        .arg("-o")
        .arg("reads")
        .args(["--min-count", "5", "-v", "-k", "9"])
        .assert()
        .success();

    let fastq_align_out_quality = Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("align")
        .arg("reads.skf")
        .output()
        .unwrap()
        .stdout;

    assert_eq!(var_hash(&fastq_align_out_quality), all_hash);
}

// Just checks that counter runs, model fit is in a unit test
#[test]
fn cov_check() {
    let sandbox = TestSetup::setup();

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("cov")
        .arg(sandbox.file_string("test_1_fwd.fastq.gz", TestDir::Input))
        .arg(sandbox.file_string("test_1_rev.fastq.gz", TestDir::Input))
        .arg("-k")
        .arg("9")
        .arg("-v")
        .assert()
        .success();

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("cov")
        .arg(sandbox.file_string("test_long_1_fwd.fastq.gz", TestDir::Input))
        .arg(sandbox.file_string("test_long_1_rev.fastq.gz", TestDir::Input))
        .arg("-k")
        .arg("33")
        .arg("-v")
        .assert()
        .success();

    // Doesn't run with fasta
    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("cov")
        .arg(sandbox.file_string("test_1.fa", TestDir::Input))
        .arg(sandbox.file_string("test_2.fa", TestDir::Input))
        .arg("-k")
        .arg("9")
        .arg("-v")
        .assert()
        .failure();
}

#[test]
fn build_auto_check() {
    let sandbox = TestSetup::setup();

    let rfile_name = sandbox.create_rfile("test", FxType::Fastq);
    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("build")
        .arg("-f")
        .arg(rfile_name)
        .arg("-o")
        .arg("reads")
        .args(["--min-count", "auto", "-v", "-k", "9", "--min-qual", "2"])
        .assert()
        .success();
}
