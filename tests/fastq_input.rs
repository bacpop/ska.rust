use snapbox::cmd::{cargo_bin, Command};

pub mod common;
use crate::common::{var_hash, TestDir, TestSetup};

// Uses perfect reads with the same fasta input as merge.skf
#[test]
fn align_fastq() {
    let sandbox = TestSetup::setup();
    let rfile_name = sandbox.create_rfile(&"test", true);

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("build")
        .arg("-f")
        .arg(rfile_name)
        .arg("-o")
        .arg("reads")
        .args(&["--min-count", "2", "-v", "-k", "9", "--min-qual", "2"])
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
        .args(&["-k", "9"])
        .arg(sandbox.file_string("test_1.fa", TestDir::Input))
        .arg(sandbox.file_string("test_2.fa", TestDir::Input))
        .args(&["-o", "fasta_k9"])
        .output()
        .unwrap()
        .stdout;

    let fasta_align_out = Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("align")
        .arg("fasta_k9.skf")
        .output()
        .unwrap()
        .stdout;

    assert_eq!(var_hash(&fastq_align_out), var_hash(&fasta_align_out));
}

// Add errors and low quality scores to split k-mers around
// CACT    TTAA    C,T
// zgrep 'TTAA.AGTG' test_1_fwd.fastq.gz
// zgrep 'CACT.TTAA' test_1_rev.fastq.gz
#[test]
fn error_fastq() {
    let sandbox = TestSetup::setup();

    // Without errors
    let rfile_name = sandbox.create_rfile(&"test", true);
    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("build")
        .arg("-f")
        .arg(rfile_name)
        .arg("-o")
        .arg("reads")
        .args(&["--min-count", "3", "-v", "-k", "9", "--min-qual", "2"])
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
    all_hash.remove(&vec!['C', 'T']);

    // With errors
    let rfile_name = sandbox.create_rfile(&"test_error", true);
    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("build")
        .arg("-f")
        .arg(rfile_name)
        .arg("-o")
        .arg("reads")
        .args(&["--min-count", "5", "-v", "-k", "9", "--min-qual", "2"])
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
    let rfile_name = sandbox.create_rfile(&"test_quality", true);
    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("build")
        .arg("-f")
        .arg(rfile_name)
        .arg("-o")
        .arg("reads")
        .args(&["--min-count", "5", "-v", "-k", "9", "--min-qual", "30"])
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
