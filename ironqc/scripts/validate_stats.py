#!/usr/bin/env python3
"""Validate samtools-style .stats SN section parity."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Dict, List, Sequence, Tuple


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
    return denom > 0.0 and (diff / denom) <= rel_tol


def find_stats_files(root: Path) -> Dict[str, Path]:
    files = {}
    for path in sorted(root.rglob("*.stats")):
        rel = str(path.relative_to(root))
        files[rel] = path
    return files


def parse_sn_section(path: Path) -> Tuple[Dict[str, str], List[str]]:
    values: Dict[str, str] = {}
    errors: List[str] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        for line_number, raw in enumerate(handle, start=1):
            if not raw.startswith("SN\t"):
                continue
            line = raw.rstrip("\r\n")
            parts = line.split("\t")
            if len(parts) != 3:
                errors.append(
                    f"{path}: line {line_number}: expected 3 tab-separated fields in SN line"
                )
                continue
            key = parts[1].strip()
            value = parts[2].strip()
            if key in values:
                errors.append(f"{path}: line {line_number}: duplicate SN key '{key}'")
            values[key] = value
    return values, errors


def compare_values(
    key: str,
    out_value: str,
    ref_value: str,
    abs_tol: float,
    rel_tol: float,
) -> str | None:
    if out_value == ref_value:
        return None

    try:
        out_num = float(out_value)
        ref_num = float(ref_value)
    except ValueError:
        return f"key '{key}': expected '{ref_value}', got '{out_value}'"

    if compare_float(out_num, ref_num, abs_tol=abs_tol, rel_tol=rel_tol):
        return None
    return (
        f"key '{key}': expected {ref_value}, got {out_value} "
        f"(abs diff={abs(out_num - ref_num):.3e})"
    )


def validate_file(
    output_path: Path, reference_path: Path, abs_tol: float, rel_tol: float
) -> List[str]:
    output_sn, output_errors = parse_sn_section(output_path)
    reference_sn, reference_errors = parse_sn_section(reference_path)
    failures: List[str] = []
    failures.extend(output_errors)
    failures.extend(reference_errors)

    out_keys = set(output_sn)
    ref_keys = set(reference_sn)

    missing = sorted(ref_keys - out_keys)
    extra = sorted(out_keys - ref_keys)
    if missing:
        failures.append(f"missing SN keys: {', '.join(missing)}")
    if extra:
        failures.append(f"extra SN keys: {', '.join(extra)}")

    for key in sorted(ref_keys & out_keys):
        maybe_error = compare_values(
            key, output_sn[key], reference_sn[key], abs_tol=abs_tol, rel_tol=rel_tol
        )
        if maybe_error:
            failures.append(maybe_error)

    return failures


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", type=Path, help="Directory with ironqc output")
    parser.add_argument(
        "reference_dir", type=Path, help="Directory with upstream reference output"
    )
    parser.add_argument(
        "--abs-tol",
        type=float,
        default=ABS_TOL,
        help=f"Absolute float tolerance (default: {ABS_TOL})",
    )
    parser.add_argument(
        "--rel-tol",
        type=float,
        default=REL_TOL,
        help=f"Relative float tolerance (default: {REL_TOL})",
    )
    args = parser.parse_args(argv)

    output_files = find_stats_files(args.output_dir)
    reference_files = find_stats_files(args.reference_dir)

    if not reference_files:
        print("FAIL: no .stats files found in reference directory")
        return 1

    failures = 0

    missing_files = sorted(set(reference_files) - set(output_files))
    extra_files = sorted(set(output_files) - set(reference_files))
    if missing_files:
        failures += len(missing_files)
        print(f"FAIL: missing output .stats files: {', '.join(missing_files)}")
    if extra_files:
        failures += len(extra_files)
        print(f"FAIL: extra output .stats files: {', '.join(extra_files)}")

    for rel_path in sorted(set(reference_files) & set(output_files)):
        file_failures = validate_file(
            output_files[rel_path],
            reference_files[rel_path],
            abs_tol=args.abs_tol,
            rel_tol=args.rel_tol,
        )
        if file_failures:
            failures += 1
            print(f"FAIL: {rel_path}")
            for item in file_failures:
                print(f"  - {item}")
        else:
            print(f"PASS: {rel_path}")

    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
