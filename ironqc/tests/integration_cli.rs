//! Integration tests for Phase 1 CLI scaffold output plumbing.

use std::fs;
use std::path::Path;
use std::process::Command;
use tempfile::tempdir;

fn touch(path: &Path) {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent).expect("failed to create parent directory");
    }
    fs::write(path, []).expect("failed to write placeholder file");
}

fn run_ironqc(args: &[&str], cwd: &Path) {
    let status = Command::new(env!("CARGO_BIN_EXE_ironqc"))
        .args(args)
        .current_dir(cwd)
        .status()
        .expect("failed to execute ironqc");

    assert!(status.success(), "ironqc command failed");
}

#[test]
fn stats_creates_expected_output() {
    let tmp = tempdir().expect("failed to create tempdir");
    let bam = tmp.path().join("sample.cram");
    let fasta = tmp.path().join("ref.fa");
    touch(&bam);
    touch(&fasta);

    run_ironqc(
        &[
            "stats",
            "sample.cram",
            "--reference",
            "ref.fa",
            "--threads",
            "2",
            "--prefix",
            "sample.md.cram",
        ],
        tmp.path(),
    );

    assert!(tmp.path().join("sample.md.cram.stats").exists());
}

#[test]
fn mosdepth_creates_expected_outputs() {
    let tmp = tempdir().expect("failed to create tempdir");
    let bam = tmp.path().join("sample.cram");
    let fasta = tmp.path().join("ref.fa");
    touch(&bam);
    touch(&fasta);

    run_ironqc(
        &[
            "mosdepth",
            "sample.cram",
            "--fasta",
            "ref.fa",
            "--threads",
            "2",
            "--by",
            "500",
            "-n",
            "--fast-mode",
            "sample.md",
        ],
        tmp.path(),
    );

    assert!(tmp
        .path()
        .join("sample.md.mosdepth.global.dist.txt")
        .exists());
    assert!(tmp
        .path()
        .join("sample.md.mosdepth.region.dist.txt")
        .exists());
    assert!(tmp.path().join("sample.md.mosdepth.summary.txt").exists());
}

#[test]
fn indexcov_creates_expected_outputs() {
    let tmp = tempdir().expect("failed to create tempdir");
    let bam_a = tmp.path().join("sample-a.bam");
    let bam_b = tmp.path().join("sample-b.bam");
    let fai = tmp.path().join("ref.fa.fai");
    touch(&bam_a);
    touch(&bam_b);
    touch(&fai);

    run_ironqc(
        &[
            "indexcov",
            "sample-a.bam",
            "sample-b.bam",
            "--fai",
            "ref.fa.fai",
            "--directory",
            "indexcov",
            "--prefix",
            "indexcov",
        ],
        tmp.path(),
    );

    assert!(tmp.path().join("indexcov/indexcov-indexcov.ped").exists());
    assert!(tmp.path().join("indexcov/indexcov-indexcov.roc").exists());
    assert!(tmp
        .path()
        .join("indexcov/indexcov-indexcov.bed.gz")
        .exists());
    assert!(tmp.path().join("indexcov/indexcov-indexcov.html").exists());
    assert!(tmp
        .path()
        .join("indexcov/indexcov-indexcov.sex.png")
        .exists());
    assert!(tmp
        .path()
        .join("indexcov/indexcov-indexcov.roc.png")
        .exists());
}

#[test]
fn bundle_creates_union_of_outputs() {
    let tmp = tempdir().expect("failed to create tempdir");
    let bam = tmp.path().join("sample.cram");
    let fasta = tmp.path().join("ref.fa");
    let fai = tmp.path().join("ref.fa.fai");
    touch(&bam);
    touch(&fasta);
    touch(&fai);

    run_ironqc(
        &[
            "bundle",
            "sample.cram",
            "--reference",
            "ref.fa",
            "--fai",
            "ref.fa.fai",
            "--threads",
            "2",
            "--prefix",
            "sample.md.cram",
            "--indexcov-dir",
            "indexcov",
        ],
        tmp.path(),
    );

    assert!(tmp.path().join("sample.md.cram.stats").exists());
    assert!(tmp
        .path()
        .join("sample.md.cram.mosdepth.global.dist.txt")
        .exists());
    assert!(tmp
        .path()
        .join("sample.md.cram.mosdepth.region.dist.txt")
        .exists());
    assert!(tmp
        .path()
        .join("sample.md.cram.mosdepth.summary.txt")
        .exists());
    assert!(tmp
        .path()
        .join("indexcov/sample.md.cram-indexcov.ped")
        .exists());
    assert!(tmp
        .path()
        .join("indexcov/sample.md.cram-indexcov.roc")
        .exists());
    assert!(tmp
        .path()
        .join("indexcov/sample.md.cram-indexcov.bed.gz")
        .exists());
    assert!(tmp
        .path()
        .join("indexcov/sample.md.cram-indexcov.html")
        .exists());
    assert!(tmp
        .path()
        .join("indexcov/sample.md.cram-indexcov.sex.png")
        .exists());
    assert!(tmp
        .path()
        .join("indexcov/sample.md.cram-indexcov.roc.png")
        .exists());
}
