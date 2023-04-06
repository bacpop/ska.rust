use std::{
    fs::{read_dir, File},
    io::{LineWriter, Write},
    path::{Path, PathBuf},
};

use assert_fs::{prelude::*, TempDir};
use predicates::prelude::*;

use hashbrown::HashSet;

#[cfg(test)]
use pretty_assertions::assert_eq;

// Creates correct path for input/output files
static FILE_IN: &'static str = "tests/test_files_in";
static FILE_TEST: &'static str = "tests/test_results_correct";
static SYM_IN: &'static str = "input";
static SYM_TEST: &'static str = "correct";
static RFILE_NAME: &'static str = "file_list.txt";

#[derive(Debug, PartialEq, Copy, Clone)]
pub enum TestDir {
    Input,
    Correct,
}

#[derive(Debug, PartialEq, Copy, Clone)]
pub enum FxType {
    Fasta,
    Fastq,
}

pub struct TestSetup {
    wd: TempDir,
}

impl TestSetup {
    pub fn setup() -> Self {
        let wd = assert_fs::TempDir::new().unwrap();
        wd.child(SYM_IN)
            .symlink_to_dir(
                Path::new(FILE_IN)
                    .canonicalize()
                    .expect("Could not link expected files"),
            )
            .unwrap();
        wd.child(SYM_TEST)
            .symlink_to_dir(
                Path::new(FILE_TEST)
                    .canonicalize()
                    .expect("Could not link expected files"),
            )
            .unwrap();
        Self { wd }
    }

    pub fn get_wd(&self) -> String {
        self.wd.path().display().to_string()
    }

    pub fn file_path(&self, name: &str, file_type: TestDir) -> PathBuf {
        match file_type {
            TestDir::Input => {
                PathBuf::from(&format!("{}/{}/{}", self.wd.path().display(), SYM_IN, name))
            }
            TestDir::Correct => PathBuf::from(&format!(
                "{}/{}/{}",
                self.wd.path().display(),
                SYM_TEST,
                name
            )),
        }
    }

    pub fn file_string(&self, name: &str, file_type: TestDir) -> String {
        self.file_path(name, file_type)
            .to_str()
            .expect("Could not unpack file path")
            .to_owned()
    }

    pub fn file_check(&self, name_out: &str, name_correct: &str) -> bool {
        let predicate_file = predicate::path::eq_file(self.wd.child(name_out).path());
        predicate_file.eval(self.file_path(name_correct, TestDir::Correct).as_path())
    }

    pub fn file_exists(&self, name_out: &str) -> bool {
        let predicate_fn = predicate::path::is_file();
        predicate_fn.eval(self.wd.child(name_out).path())
    }

    pub fn create_rfile(&self, prefix: &str, fx_type: FxType) -> &str {
        // Create an rfile in the tmp dir
        let mut rfile = LineWriter::new(
            File::create(format!("{}/{}", self.get_wd(), RFILE_NAME))
                .expect("Could not write rfile"),
        );
        match fx_type {
            FxType::Fastq => {
                writeln!(
                    rfile,
                    "{}",
                    &format!(
                        "{}_1\t{}\t{}",
                        prefix,
                        self.file_string(&format!("{}_1_fwd.fastq.gz", prefix), TestDir::Input),
                        self.file_string(&format!("{}_1_rev.fastq.gz", prefix), TestDir::Input),
                    )
                )
                .unwrap();
                writeln!(
                    rfile,
                    "{}",
                    &format!(
                        "{}_2\t{}\t{}",
                        prefix,
                        self.file_string(&format!("{}_2_fwd.fastq.gz", prefix), TestDir::Input),
                        self.file_string(&format!("{}_2_rev.fastq.gz", prefix), TestDir::Input),
                    )
                )
                .unwrap();
            }
            FxType::Fasta => {
                writeln!(
                    rfile,
                    "{}",
                    &format!(
                        "{}_1\t{}",
                        prefix,
                        self.file_string(&format!("{}_1.fa", prefix), TestDir::Input),
                    )
                )
                .unwrap();
                writeln!(
                    rfile,
                    "{}",
                    &format!(
                        "{}_2\t{}",
                        prefix,
                        self.file_string(&format!("{}_2.fa", prefix), TestDir::Input),
                    )
                )
                .unwrap();
            }
        };
        RFILE_NAME
    }

    pub fn create_par_rfile(&self) -> &str {
        let mut rfile = LineWriter::new(
            File::create(format!("{}/{}", self.get_wd(), RFILE_NAME))
                .expect("Could not write rfile"),
        );

        let paths = read_dir(self.file_string("par_test", TestDir::Input)).unwrap();
        for path in paths {
            let path_str = path.unwrap().path().display().to_string();
            writeln!(rfile, "{path_str}\t{path_str}").unwrap();
        }
        RFILE_NAME
    }
}

// Helper for multi sample fasta, where column order isn't conserved
pub fn var_hash(aln_string: &[u8]) -> HashSet<Vec<char>> {
    let fastq_align_out = String::from_utf8(aln_string.to_vec()).unwrap();
    let mut sample_vecs: Vec<Vec<char>> = Vec::new();

    // Read sample lines
    for (idx, line) in fastq_align_out.lines().enumerate() {
        if idx % 2 == 1 {
            let line_vec: Vec<char> = line.chars().collect();
            sample_vecs.push(line_vec);
        }
    }

    // Transpose into a set
    let mut variant_pairs: HashSet<Vec<char>> = HashSet::new();
    for var_idx in 0..sample_vecs[0].len() {
        let mut var_vec = Vec::new();
        for sample_vec in &sample_vecs {
            var_vec.push(sample_vec[var_idx]);
        }
        variant_pairs.insert(var_vec);
    }

    return variant_pairs;
}

// Helper for comparing mapped alignments with different sample names
pub fn cmp_map_aln(aln1: &[u8], aln2: &[u8]) {
    let aln1_str = String::from_utf8(aln1.to_vec()).unwrap();
    let aln2_str = String::from_utf8(aln2.to_vec()).unwrap();

    for (line1, line2) in aln1_str.lines().zip(aln2_str.lines()).skip(1).step_by(2) {
        assert_eq!(line1, line2);
    }
}

// Helper for checking alignment length
pub fn aln_length(aln1: &[u8]) -> Vec<usize> {
    let aln1_str = String::from_utf8(aln1.to_vec()).unwrap();

    let mut lengths = Vec::new();
    for aln_line in aln1_str.lines().skip(1).step_by(2) {
        lengths.push(aln_line.len());
    }
    lengths
}
