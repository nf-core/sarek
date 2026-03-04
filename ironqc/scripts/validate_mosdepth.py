#!/usr/bin/env python3
"""Validate mosdepth summary/global/region outputs against reference."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Dict, List, Sequence, Tuple


ABS_TOL = 1e-10
REL_TOL = 1e-6

SUFFIXES = (
    ".mosdepth.summary.txt",
    ".mosdepth.global.dist.txt",
    ".mosdepth.region.dist.txt",
)


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


def discover_files(root: Path, suffix: str) -> Dict[str, Path]:
    return {
        str(path.relative_to(root)): path for path in sorted(root.rglob(f"*{suffix}"))
    }


def parse_tsv(path: Path) -> List[List[str]]:
    rows: List[List[str]] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.rstrip("\r\n")
            if not line:
                rows.append([])
            else:
                rows.append(line.split("\t"))
            if raw.endswith(" "):
                raise ValueError(f"{path}: line {line_number}: trailing whitespace")
    return rows


def _compare_dist_rows(
    output_rows: List[List[str]],
    reference_rows: List[List[str]],
    abs_tol: float,
    rel_tol: float,
) -> List[str]:
    failures: List[str] = []
    if len(output_rows) != len(reference_rows):
        failures.append(
            f"row count mismatch: expected {len(reference_rows)}, got {len(output_rows)}"
        )
        return failures

    for idx, (out_row, ref_row) in enumerate(zip(output_rows, reference_rows), start=1):
        if len(out_row) != 3 or len(ref_row) != 3:
            failures.append(f"line {idx}: expected 3 columns for dist rows")
            continue
        if out_row[0] != ref_row[0]:
            failures.append(
                f"line {idx}, col 1: expected '{ref_row[0]}', got '{out_row[0]}'"
            )
        if out_row[1] != ref_row[1]:
            failures.append(
                f"line {idx}, col 2: expected '{ref_row[1]}', got '{out_row[1]}'"
            )

        try:
            out_val = float(out_row[2])
            ref_val = float(ref_row[2])
        except ValueError:
            failures.append(f"line {idx}, col 3: non-numeric value")
            continue

        if not compare_float(out_val, ref_val, abs_tol, rel_tol):
            failures.append(
                f"line {idx}, col 3: expected {ref_row[2]}, got {out_row[2]} "
                f"(abs diff={abs(out_val - ref_val):.3e})"
            )
    return failures


def _compare_summary_rows(
    output_rows: List[List[str]],
    reference_rows: List[List[str]],
    abs_tol: float,
    rel_tol: float,
) -> List[str]:
    failures: List[str] = []
    if len(output_rows) != len(reference_rows):
        failures.append(
            f"row count mismatch: expected {len(reference_rows)}, got {len(output_rows)}"
        )
        return failures
    if not reference_rows:
        failures.append("empty reference summary")
        return failures

    if output_rows[0] != reference_rows[0]:
        failures.append(
            "header mismatch: "
            f"expected {'|'.join(reference_rows[0])}, got {'|'.join(output_rows[0])}"
        )

    for idx, (out_row, ref_row) in enumerate(
        zip(output_rows[1:], reference_rows[1:]), start=2
    ):
        if len(out_row) != 6 or len(ref_row) != 6:
            failures.append(f"line {idx}: expected 6 columns for summary rows")
            continue
        for col_idx in (0, 1, 2, 4, 5):
            if out_row[col_idx] != ref_row[col_idx]:
                failures.append(
                    f"line {idx}, col {col_idx + 1}: expected '{ref_row[col_idx]}', got '{out_row[col_idx]}'"
                )

        try:
            out_mean = float(out_row[3])
            ref_mean = float(ref_row[3])
        except ValueError:
            failures.append(f"line {idx}, col 4: non-numeric mean value")
            continue
        if not compare_float(out_mean, ref_mean, abs_tol, rel_tol):
            failures.append(
                f"line {idx}, col 4: expected {ref_row[3]}, got {out_row[3]} "
                f"(abs diff={abs(out_mean - ref_mean):.3e})"
            )
    return failures


def compare_file(
    output_path: Path, reference_path: Path, abs_tol: float, rel_tol: float
) -> List[str]:
    output_rows = parse_tsv(output_path)
    reference_rows = parse_tsv(reference_path)
    if output_path.name.endswith("summary.txt"):
        return _compare_summary_rows(output_rows, reference_rows, abs_tol, rel_tol)
    return _compare_dist_rows(output_rows, reference_rows, abs_tol, rel_tol)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("reference_dir", type=Path)
    parser.add_argument("--abs-tol", type=float, default=ABS_TOL)
    parser.add_argument("--rel-tol", type=float, default=REL_TOL)
    args = parser.parse_args(argv)

    total_failures = 0

    for suffix in SUFFIXES:
        output_files = discover_files(args.output_dir, suffix)
        reference_files = discover_files(args.reference_dir, suffix)

        missing = sorted(set(reference_files) - set(output_files))
        extra = sorted(set(output_files) - set(reference_files))
        if missing:
            total_failures += len(missing)
            print(f"FAIL [{suffix}]: missing files: {', '.join(missing)}")
        if extra:
            total_failures += len(extra)
            print(f"FAIL [{suffix}]: extra files: {', '.join(extra)}")

        for rel in sorted(set(reference_files) & set(output_files)):
            failures = compare_file(
                output_files[rel],
                reference_files[rel],
                abs_tol=args.abs_tol,
                rel_tol=args.rel_tol,
            )
            if failures:
                total_failures += 1
                print(f"FAIL: {rel}")
                for item in failures:
                    print(f"  - {item}")
            else:
                print(f"PASS: {rel}")

    return 1 if total_failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
