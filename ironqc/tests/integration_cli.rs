use std::fs;
use std::path::Path;
use std::process::Command;
use tempfile::tempdir;

fn create_file(path: &Path, content: &[u8]) {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent).expect("create parent dir");
    }
    fs::write(path, content).expect("write file");
}

fn run_ironqc(args: &[&str], cwd: &Path) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_ironqc"))
        .args(args)
        .current_dir(cwd)
        .output()
        .expect("failed to execute ironqc")
}

#[test]
fn indexcov_creates_expected_outputs() {
    let tmp = tempdir().expect("tempdir");
    let bam_a = tmp.path().join("sample-a.bam");
    let bam_b = tmp.path().join("sample-b.bam");
    let fai = tmp.path().join("ref.fa.fai");

    create_file(&bam_a, b"");
    create_file(&bam_b, b"");
    create_file(
        &fai,
        b"chr1\t1000000\t6\t80\t81\nchr2\t500000\t1000100\t80\t81\n",
    );

    let output = run_ironqc(
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

    assert!(
        output.status.success(),
        "indexcov failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(tmp.path().join("indexcov/indexcov-indexcov.ped").exists());
    assert!(tmp.path().join("indexcov/indexcov-indexcov.roc").exists());
    assert!(tmp
        .path()
        .join("indexcov/indexcov-indexcov.bed.gz")
        .exists());
    assert!(tmp.path().join("indexcov/indexcov-indexcov.html").exists());
}

#[test]
fn indexcov_ped_has_correct_header_and_samples() {
    let tmp = tempdir().expect("tempdir");
    let bam = tmp.path().join("mysample.bam");
    let fai = tmp.path().join("ref.fa.fai");

    create_file(&bam, b"");
    create_file(&fai, b"chr1\t248956422\t6\t80\t81\n");

    let output = run_ironqc(
        &[
            "indexcov",
            "mysample.bam",
            "--fai",
            "ref.fa.fai",
            "--directory",
            "out",
            "--prefix",
            "test",
        ],
        tmp.path(),
    );
    assert!(output.status.success());

    let ped_content =
        fs::read_to_string(tmp.path().join("out/test-indexcov.ped")).expect("read ped");
    assert!(ped_content.starts_with("#family_id\tsample_id"));
    assert!(ped_content.contains("mysample"));
    assert!(ped_content.contains("CN_chr1"));
}

#[test]
fn indexcov_roc_has_correct_structure() {
    let tmp = tempdir().expect("tempdir");
    let bam = tmp.path().join("s1.bam");
    let fai = tmp.path().join("ref.fa.fai");

    create_file(&bam, b"");
    create_file(&fai, b"chr1\t1000000\t6\t80\t81\n");

    let output = run_ironqc(
        &[
            "indexcov",
            "s1.bam",
            "--fai",
            "ref.fa.fai",
            "--directory",
            "roc_test",
            "--prefix",
            "roc",
        ],
        tmp.path(),
    );
    assert!(output.status.success());

    let roc_content =
        fs::read_to_string(tmp.path().join("roc_test/roc-indexcov.roc")).expect("read roc");
    assert!(roc_content.starts_with("#chrom\tcov"));
    assert!(roc_content.contains("chr1"));
    assert!(roc_content.contains("s1"));
}

#[test]
fn stats_fails_gracefully_on_invalid_bam() {
    let tmp = tempdir().expect("tempdir");
    let bam = tmp.path().join("bad.bam");
    let fasta = tmp.path().join("ref.fa");

    create_file(&bam, b"not a bam file");
    create_file(&fasta, b">chr1\nACGT\n");

    let output = run_ironqc(
        &[
            "stats",
            "bad.bam",
            "--reference",
            "ref.fa",
            "--prefix",
            "test",
        ],
        tmp.path(),
    );

    assert!(!output.status.success());
}

#[test]
fn mosdepth_fails_gracefully_on_invalid_bam() {
    let tmp = tempdir().expect("tempdir");
    let bam = tmp.path().join("bad.bam");
    let fasta = tmp.path().join("ref.fa");

    create_file(&bam, b"not a bam file");
    create_file(&fasta, b">chr1\nACGT\n");

    let output = run_ironqc(
        &[
            "mosdepth", "bad.bam", "--fasta", "ref.fa", "--by", "500", "test",
        ],
        tmp.path(),
    );

    assert!(!output.status.success());
}
