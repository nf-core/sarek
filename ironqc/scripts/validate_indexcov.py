#!/usr/bin/env python3
"""Validate indexcov PED and ROC output parity."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Dict, List, Sequence


ABS_TOL = 1e-10
REL_TOL = 1e-6


def compare_float(a: float, b: float, abs_tol: float, rel_tol: float) -> bool:
    if a == b:
        return True
    if math.isnan(a) and math.isnan(b):
        return True
    diff = abs(a - b)
    if diff <= abs_tol:
        return True
    denom = max(abs(a), abs(b))
    return denom > 0 and (diff / denom) <= rel_tol


def discover(root: Path, suffix: str) -> Dict[str, Path]:
    return {
        str(path.relative_to(root)): path for path in sorted(root.rglob(f"*{suffix}"))
    }


def read_tsv(path: Path) -> List[List[str]]:
    rows: List[List[str]] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.rstrip("\r\n")
            if not line:
                rows.append([])
                continue
            if line.endswith(" "):
                raise ValueError(f"{path}: line {line_number}: trailing whitespace")
            rows.append(line.split("\t"))
    return rows


def _compare_line(
    out_row: List[str],
    ref_row: List[str],
    line_number: int,
    numeric_start_col: int,
    abs_tol: float,
    rel_tol: float,
) -> List[str]:
    failures: List[str] = []
    if len(out_row) != len(ref_row):
        failures.append(
            f"line {line_number}: column count mismatch expected {len(ref_row)}, got {len(out_row)}"
        )
        return failures

    for idx, (out_val, ref_val) in enumerate(zip(out_row, ref_row), start=1):
        if idx < numeric_start_col:
            if out_val != ref_val:
                failures.append(
                    f"line {line_number}, col {idx}: expected '{ref_val}', got '{out_val}'"
                )
            continue

        if out_val == ref_val:
            continue
        if out_val == "NA" or ref_val == "NA":
            failures.append(
                f"line {line_number}, col {idx}: expected '{ref_val}', got '{out_val}'"
            )
            continue
        try:
            out_num = float(out_val)
            ref_num = float(ref_val)
        except ValueError:
            failures.append(
                f"line {line_number}, col {idx}: expected '{ref_val}', got '{out_val}'"
            )
            continue

        if not compare_float(out_num, ref_num, abs_tol, rel_tol):
            failures.append(
                f"line {line_number}, col {idx}: expected {ref_val}, got {out_val} "
                f"(abs diff={abs(out_num - ref_num):.3e})"
            )

    return failures


def compare_ped(
    output_path: Path, reference_path: Path, abs_tol: float, rel_tol: float
) -> List[str]:
    output_rows = read_tsv(output_path)
    reference_rows = read_tsv(reference_path)
    failures: List[str] = []
    if len(output_rows) != len(reference_rows):
        failures.append(
            f"row count mismatch: expected {len(reference_rows)}, got {len(output_rows)}"
        )
        return failures
    if not reference_rows:
        failures.append("empty reference PED")
        return failures

    if output_rows[0] != reference_rows[0]:
        failures.append("header mismatch")

    for line_number, (out_row, ref_row) in enumerate(
        zip(output_rows[1:], reference_rows[1:]), start=2
    ):
        failures.extend(
            _compare_line(
                out_row,
                ref_row,
                line_number,
                numeric_start_col=7,
                abs_tol=abs_tol,
                rel_tol=rel_tol,
            )
        )
    return failures


def compare_roc(
    output_path: Path, reference_path: Path, abs_tol: float, rel_tol: float
) -> List[str]:
    output_rows = read_tsv(output_path)
    reference_rows = read_tsv(reference_path)
    failures: List[str] = []
    if len(output_rows) != len(reference_rows):
        failures.append(
            f"row count mismatch: expected {len(reference_rows)}, got {len(output_rows)}"
        )
        return failures
    if not reference_rows:
        failures.append("empty reference ROC")
        return failures

    if output_rows[0] != reference_rows[0]:
        failures.append("header mismatch")

    for line_number, (out_row, ref_row) in enumerate(
        zip(output_rows[1:], reference_rows[1:]), start=2
    ):
        failures.extend(
            _compare_line(
                out_row,
                ref_row,
                line_number,
                numeric_start_col=2,
                abs_tol=abs_tol,
                rel_tol=rel_tol,
            )
        )
    return failures


def run_group(
    output_dir: Path,
    reference_dir: Path,
    suffix: str,
    comparator,
    abs_tol: float,
    rel_tol: float,
) -> int:
    failures = 0
    output_files = discover(output_dir, suffix)
    reference_files = discover(reference_dir, suffix)

    missing = sorted(set(reference_files) - set(output_files))
    extra = sorted(set(output_files) - set(reference_files))
    if missing:
        failures += len(missing)
        print(f"FAIL [{suffix}]: missing files: {', '.join(missing)}")
    if extra:
        failures += len(extra)
        print(f"FAIL [{suffix}]: extra files: {', '.join(extra)}")

    for rel_path in sorted(set(reference_files) & set(output_files)):
        file_failures = comparator(
            output_files[rel_path],
            reference_files[rel_path],
            abs_tol,
            rel_tol,
        )
        if file_failures:
            failures += 1
            print(f"FAIL: {rel_path}")
            for item in file_failures:
                print(f"  - {item}")
        else:
            print(f"PASS: {rel_path}")

    return failures


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("reference_dir", type=Path)
    parser.add_argument("--abs-tol", type=float, default=ABS_TOL)
    parser.add_argument("--rel-tol", type=float, default=REL_TOL)
    args = parser.parse_args(argv)

    failures = 0
    failures += run_group(
        args.output_dir,
        args.reference_dir,
        "-indexcov.ped",
        compare_ped,
        args.abs_tol,
        args.rel_tol,
    )
    failures += run_group(
        args.output_dir,
        args.reference_dir,
        "-indexcov.roc",
        compare_roc,
        args.abs_tol,
        args.rel_tol,
    )

    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
