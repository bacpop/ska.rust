use snapbox::cmd::{cargo_bin, Command};

#[cfg(test)]
use pretty_assertions::assert_eq;

pub mod common;
use crate::common::*;

// NB: to view output, uncomment the current_dir lines

#[test]
fn map_aln() {
    let sandbox = TestSetup::setup();

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("map")
        .arg(sandbox.file_string("test_ref.fa", TestDir::Input))
        .arg(sandbox.file_string("test_1.fa", TestDir::Input))
        .arg(sandbox.file_string("test_2.fa", TestDir::Input))
        .arg("-o")
        .arg("map.aln")
        .assert()
        .success();

    assert_eq!(true, sandbox.file_exists("map.aln"));

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("map")
        .arg(sandbox.file_string("test_ref.fa", TestDir::Input))
        .arg(sandbox.file_string("merge.skf", TestDir::Input))
        .args(&["-v", "--threads", "2"])
        .assert()
        .stdout_eq_path(sandbox.file_string("map_aln.stdout", TestDir::Correct));

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("map")
        .arg(sandbox.file_string("test_ref.fa", TestDir::Input))
        .arg(sandbox.file_string("merge_k9.skf", TestDir::Input))
        .assert()
        .stdout_eq_path(sandbox.file_string("map_aln_k9.stdout", TestDir::Correct));

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("map")
        .arg(sandbox.file_string("test_ref.fa", TestDir::Input))
        .arg(sandbox.file_string("merge_k9.skf", TestDir::Input))
        .arg("--ambig-mask")
        .arg("-v")
        .assert()
        .stdout_eq_path(sandbox.file_string("map_aln_k9_filter.stdout", TestDir::Correct));

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("map")
        .arg(sandbox.file_string("test_ref_two_chrom.fa", TestDir::Input))
        .arg(sandbox.file_string("merge.skf", TestDir::Input))
        .assert()
        .stdout_eq_path(sandbox.file_string("map_aln_two_chrom.stdout", TestDir::Correct));

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("map")
        .arg(sandbox.file_string("test_ref.fa", TestDir::Input))
        .arg(sandbox.file_string("test_1.fa", TestDir::Input))
        .arg(sandbox.file_string("indel_test.fa", TestDir::Input))
        .assert()
        .stdout_eq_path(sandbox.file_string("map_aln_indels.stdout", TestDir::Correct));

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("build")
        .arg("-k")
        .arg("17")
        .arg(sandbox.file_string("ambig_test_1.fa", TestDir::Input))
        .arg(sandbox.file_string("ambig_test_2.fa", TestDir::Input))
        .arg("-o")
        .arg("ambig_map")
        .arg("--single-strand")
        .assert()
        .success();

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("map")
        .arg(sandbox.file_string("ambig_test_ref.fa", TestDir::Input))
        .arg("ambig_map.skf")
        .assert()
        .stdout_eq_path(sandbox.file_string("map_aln_ambig.stdout", TestDir::Correct));
}

#[test]
fn map_u128() {
    let sandbox = TestSetup::setup();

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("map")
        .arg(sandbox.file_string("test_ref.fa", TestDir::Input))
        .arg(sandbox.file_string("merge_k41.skf", TestDir::Input))
        .args(&["-v", "--threads", "2"])
        .assert()
        .stdout_eq_path(sandbox.file_string("map_aln_k41.stdout", TestDir::Correct));

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("map")
        .arg(sandbox.file_string("test_ref.fa", TestDir::Input))
        .arg(sandbox.file_string("merge_k41.skf", TestDir::Input))
        .args(&["-v", "--threads", "2"])
        .arg("-f")
        .arg("vcf")
        .assert()
        .stdout_eq_path(sandbox.file_string("map_vcf_k41.stdout", TestDir::Correct));
}

#[test]
fn map_vcf() {
    let sandbox = TestSetup::setup();

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("map")
        .arg(sandbox.file_string("test_ref.fa", TestDir::Input))
        .arg(sandbox.file_string("test_1.fa", TestDir::Input))
        .arg(sandbox.file_string("test_2.fa", TestDir::Input))
        .arg("-o")
        .arg("map.vcf")
        .arg("-f")
        .arg("vcf")
        .assert()
        .success();

    assert_eq!(true, sandbox.file_exists("map.vcf"));

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("map")
        .arg(sandbox.file_string("test_ref.fa", TestDir::Input))
        .arg(sandbox.file_string("merge.skf", TestDir::Input))
        .arg("-f")
        .arg("vcf")
        .assert()
        .stdout_eq_path(sandbox.file_string("map_vcf.stdout", TestDir::Correct));

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("map")
        .arg(sandbox.file_string("test_ref_two_chrom.fa", TestDir::Input))
        .arg(sandbox.file_string("merge.skf", TestDir::Input))
        .arg("-f")
        .arg("vcf")
        .assert()
        .stdout_eq_path(sandbox.file_string("map_vcf_two_chrom.stdout", TestDir::Correct));

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("map")
        .arg(sandbox.file_string("test_ref.fa", TestDir::Input))
        .arg(sandbox.file_string("test_1.fa", TestDir::Input))
        .arg(sandbox.file_string("indel_test.fa", TestDir::Input))
        .arg("-f")
        .arg("vcf")
        .assert()
        .stdout_eq_path(sandbox.file_string("map_vcf_indels.stdout", TestDir::Correct));
}

#[test]
fn map_rev_comp() {
    let sandbox = TestSetup::setup();

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("build")
        .arg("-o")
        .arg("rc_build")
        .arg("-k")
        .arg("9")
        .arg(sandbox.file_string("test_1.fa", TestDir::Input))
        .arg(sandbox.file_string("test_2_rc.fa", TestDir::Input))
        .assert()
        .success();

    let rc_map = Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("map")
        .arg(sandbox.file_string("test_ref.fa", TestDir::Input))
        .arg("fwd_build.skf")
        .output()
        .unwrap()
        .stdout;

    let fwd_map = Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("map")
        .arg(sandbox.file_string("test_ref.fa", TestDir::Input))
        .arg(sandbox.file_string("merge.skf", TestDir::Input))
        .output()
        .unwrap()
        .stdout;

    cmp_map_aln(&rc_map, &fwd_map);

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("build")
        .arg("-o")
        .arg("ss_map")
        .arg("-k")
        .arg("9")
        .arg("--single-strand")
        .arg(sandbox.file_string("test_1.fa", TestDir::Input))
        .arg(sandbox.file_string("test_2_rc.fa", TestDir::Input))
        .assert()
        .success();

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("map")
        .arg(sandbox.file_string("test_ref.fa", TestDir::Input))
        .arg("ss_map.skf")
        .assert()
        .stdout_eq_path(sandbox.file_string("map_ss.stdout", TestDir::Correct));

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("map")
        .arg(sandbox.file_string("test_ref.fa", TestDir::Input))
        .arg("ss_map.skf")
        .arg("-f")
        .arg("vcf")
        .assert()
        .stdout_eq_path(sandbox.file_string("map_vcf_ss.stdout", TestDir::Correct));
}

// Tests the --repeat-mask option
#[test]
fn repeat_mask() {
    let sandbox = TestSetup::setup();

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("map")
        .arg(sandbox.file_string("test_ref.fa", TestDir::Input))
        .arg(sandbox.file_string("merge_k9.skf", TestDir::Input))
        .arg("-v")
        .arg("--repeat-mask")
        .assert()
        .stdout_eq_path(sandbox.file_string("map_aln_k9.masked.stdout", TestDir::Correct));

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("map")
        .arg(sandbox.file_string("test_ref.fa", TestDir::Input))
        .arg(sandbox.file_string("merge_k9.skf", TestDir::Input))
        .arg("-v")
        .arg("--repeat-mask")
        .args(&["--format", "vcf"])
        .assert()
        .stdout_eq_path(sandbox.file_string("map_vcf_k9.masked.stdout", TestDir::Correct));

    // Two identical chromosomes. All bases masked
    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("map")
        .arg(sandbox.file_string("test_ref_two_chrom.fa", TestDir::Input))
        .arg(sandbox.file_string("merge_k9.skf", TestDir::Input))
        .arg("-v")
        .arg("--repeat-mask")
        .assert()
        .stdout_eq_path(sandbox.file_string("map_all_repeats.masked.stdout", TestDir::Correct));

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("map")
        .arg(sandbox.file_string("test_ref_two_chrom_repeats.fa", TestDir::Input))
        .arg(sandbox.file_string("merge_k9.skf", TestDir::Input))
        .arg("-v")
        .arg("--repeat-mask")
        .assert()
        .stdout_eq_path(sandbox.file_string("map_aln_two_chrom.masked.stdout", TestDir::Correct));

    Command::new(cargo_bin("ska"))
        .current_dir(sandbox.get_wd())
        .arg("map")
        .arg(sandbox.file_string("test_ref_two_chrom_repeats.fa", TestDir::Input))
        .arg(sandbox.file_string("merge_k9.skf", TestDir::Input))
        .arg("-v")
        .arg("--repeat-mask")
        .args(&["--format", "vcf"])
        .assert()
        .stdout_eq_path(sandbox.file_string("map_vcf_two_chrom.masked.stdout", TestDir::Correct));
}
