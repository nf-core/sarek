from pathlib import Path

import pytest

from checks.if_to_filter import check_source, main
from checks.nextflow import nextflow_files

# No warning expected. This does not mean the example is preferred style.
VALID = {
    "already filtered": """
workflow {
    ALIGN(reads.filter { _read -> params.align }, reference)
}
""",
    "else block": """
workflow {
    if (params.align) {
        ALIGN(reads)
    } else {
        OTHER_ALIGN(reads)
    }
}
""",
    "else if": """
workflow {
    if (params.align) { ALIGN(reads) }
    else if (params.other) { OTHER_ALIGN(reads) }
}
""",
    "call plus assignment": """
workflow {
    if (params.align) {
        ALIGN(reads)
        output = ALIGN.out.bam
    }
}
""",
    "two calls": """
workflow {
    if (params.align) {
        ALIGN(reads)
        INDEX(ALIGN.out.bam)
    }
}
""",
    "assignment then call": """
workflow {
    if (params.align) {
        reads = channel.of('reads.fastq')
        ALIGN(reads)
    }
}
""",
    "unbraced guard": """
workflow {
    if (params.align) ALIGN(reads)
}
""",
    "lowercase call": """
workflow {
    if (!params.input) { error('Missing input') }
}
""",
    "inside closure": """
workflow {
    reads.map { read ->
        if (params.align) { ALIGN(read) }
    }
}
""",
    "inside function": """
def choose(reads) {
    if (params.align) { ALIGN(reads) }
}
""",
    "qualified call": """
workflow {
    if (params.align) { helpers.ALIGN(reads) }
}
""",
    "assigned call": """
workflow {
    if (params.align) { output = ALIGN(reads) }
}
""",
    "nested guard": """
workflow {
    if (params.outer) {
        if (params.align) { ALIGN(reads) }
    }
}
""",
    "suppressed": """
workflow {
    // ast-grep-ignore: if-to-filter
    if (params.align) { ALIGN(reads) }
}
""",
}

# (source, process, 1-based line, 1-based column of the reported `if`)
INVALID = {
    "inline block": (
        """workflow {
    if (params.align) { ALIGN(reads) }
}
""",
        "ALIGN",
        2,
        5,
    ),
    "main section": (
        """workflow ALIGNMENT {
    take:
    reads
    reference
    main:
    if (params.align) {
        ALIGN(reads, reference)
    }
}
""",
        "ALIGN",
        6,
        5,
    ),
    "complex condition and closure argument": (
        """workflow {
    if (!params.skip_alignment && params.tools?.split(',')?.contains('align')) {
        ALIGN_FASTQ(reads.map { meta, fastq -> [meta, fastq] }, reference)
    }
}
""",
        "ALIGN_FASTQ",
        2,
        5,
    ),
    "comments around call": (
        """workflow ALIGNMENT {
    // Run the selected aligner.
    if (params.align) {
        // Keep the reference input reusable.
        ALIGN(reads, reference)
        /* Reference remains a value channel. */
    }
}
""",
        "ALIGN",
        3,
        5,
    ),
    "suppression for another rule": (
        """workflow {
    // ast-grep-ignore: some-other-rule
    if (params.align) { ALIGN(reads) }
}
""",
        "ALIGN",
        3,
        5,
    ),
}


@pytest.mark.parametrize("source", VALID.values(), ids=VALID.keys())
def test_valid(source: str) -> None:
    assert list(check_source(source)) == []


@pytest.mark.parametrize(
    ("source", "process", "line", "column"), INVALID.values(), ids=INVALID.keys()
)
def test_invalid(source: str, process: str, line: int, column: int) -> None:
    (finding,) = check_source(source)
    start, end = finding.node.range().start, finding.node.range().end
    assert (start.line + 1, start.column + 1) == (line, column)
    assert finding.node.kind() == "if_statement"
    assert finding.node.text().startswith("if (") and finding.node.text().endswith("}")
    assert source[start.index : end.index] == finding.node.text()
    assert f"conditionally invoking {process}." in finding.message


def test_ignored_directories(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    root = tmp_path / "work" / "sarek"  # an ignored name above the scan root must not matter
    source = INVALID["inline block"][0]
    kept = ["main.nf", "workflows/sarek/main.nf", "subworkflows/local/x/main.nf"]
    ignored = [
        "modules/local/x/main.nf",
        "subworkflows/nf-core/x/main.nf",
        "tests/main.nf",
        "work/ab/main.nf",
        ".nf-test/main.nf",
    ]
    for name in kept + ignored:
        (root / name).parent.mkdir(parents=True, exist_ok=True)
        (root / name).write_text(source)
    monkeypatch.chdir(root)
    found = {p.relative_to(root).as_posix() for p in nextflow_files([root])}
    assert found == set(kept)


def test_cli_exit_codes(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    script = tmp_path / "main.nf"
    script.write_text(INVALID["inline block"][0])
    assert main([str(script)]) == 0
    assert main(["--error", str(script)]) == 1
    assert f"{script}:2:5" in capsys.readouterr().out

    script.write_text(VALID["two calls"])
    assert main(["--error", str(script)]) == 0
