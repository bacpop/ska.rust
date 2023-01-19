
use std::{
    fs::File,
    io::{LineWriter, Write},
};

use snapbox::cmd::{cargo_bin, Command};

mod common;
use crate::common::{TestDir, TestSetup};

use hashbrown::HashSet;

// Helper for two sample fasta
fn var_hash(aln_string: &[u8]) -> HashSet<(char, char)> {
    let fastq_align_out = String::from_utf8(aln_string.to_vec()).unwrap();
    let mut sample1: Vec<char> = Vec::new();
    let mut variant_pairs: HashSet<(char, char)> = HashSet::new();
    for (idx, line) in fastq_align_out.lines().enumerate() {
        if idx == 1 {
            sample1 = line.chars().collect();
        } else if idx == 3 {
            for (first, second) in sample1.iter().zip(line.chars()) {
                variant_pairs.insert((*first, second));
            }
        }
    }
    return variant_pairs;
}

// Uses perfect reads with the same fasta input as merge.skf
#[test]
fn align_fastq() {
    let sandbox = TestSetup::setup();

    // Create an rfile in the tmp dir
    let rfile_name = "file_list.txt";
    let mut rfile = LineWriter::new(
        File::create(format!("{}/{}", sandbox.get_wd(), rfile_name))
            .expect("Could not write rfile"),
    );
    writeln!(
        rfile,
        "{}",
        &format!(
            "test_1\t{}\t{}",
            sandbox.file_string("test_1_fwd.fastq.gz", TestDir::Input),
            sandbox.file_string("test_1_rev.fastq.gz", TestDir::Input)
        )
    )
    .unwrap();
    writeln!(
        rfile,
        "{}",
        &format!(
            "test_2\t{}\t{}",
            sandbox.file_string("test_2_fwd.fastq.gz", TestDir::Input),
            sandbox.file_string("test_2_rev.fastq.gz", TestDir::Input)
        )
    )
    .unwrap();

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("build")
        .arg("-f")
        .arg("file_list.txt")
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


