use snapbox::cmd::{cargo_bin, Command};

mod common;
use crate::common::{TestDir, TestSetup, var_hash};

// Uses perfect reads with the same fasta input as merge.skf
#[test]
fn align_fastq() {
    let sandbox = TestSetup::setup();
    let rfile_name = sandbox.create_rfile(true);

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
        .unwrap().stdout;

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("build")
        .args(&["-k", "9"])
        .arg(sandbox.file_string("test_1.fa", TestDir::Input))
        .arg(sandbox.file_string("test_2.fa", TestDir::Input))
        .args(&["-o", "fasta_k9"])
        .output()
        .unwrap().stdout;

    let fasta_align_out = Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("align")
        .arg("fasta_k9.skf")
        .output()
        .unwrap().stdout;

    assert_eq!(var_hash(&fastq_align_out), var_hash(&fasta_align_out));
}


