use assert_fs::{prelude::*, TempDir};
use predicates::prelude::*;
use std::path::{Path, PathBuf};

// Creates correct path for input/output files
static FILE_IN: &'static str = "tests/test_files_in";
static FILE_TEST: &'static str = "tests/test_results_correct";
static SYM_IN: &'static str = "input";
static SYM_TEST: &'static str = "correct";

#[derive(Debug, PartialEq, Copy, Clone)]
pub enum TestDir {
    Input,
    Correct,
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
}
