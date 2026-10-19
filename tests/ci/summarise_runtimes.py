#!/usr/bin/env python3

import argparse
import csv
import statistics
from collections import defaultdict
from pathlib import Path


def read_reports(report_dir):
    observations = defaultdict(list)
    for report in sorted(Path(report_dir).rglob("*.csv")):
        with report.open(newline="") as handle:
            for row in csv.DictReader(handle):
                if row.get("result") != "PASSED":
                    continue
                filename = Path(row["filename"])
                try:
                    filename = filename.relative_to(Path.cwd())
                except ValueError:
                    parts = filename.parts
                    if "tests" in parts:
                        filename = Path(*parts[parts.index("tests") :])
                key = f"{filename}::{row['test']}"
                observations[key].append(float(row["time"]))
    return observations


def summarise(observations):
    rows = []
    for key, durations in observations.items():
        ordered = sorted(durations)
        p75_index = max(0, int(0.75 * len(ordered) + 0.999999) - 1)
        rows.append(
            {
                "test": key,
                "observations": len(ordered),
                "median_seconds": round(statistics.median(ordered), 3),
                "p75_seconds": round(ordered[p75_index], 3),
                "max_seconds": round(max(ordered), 3),
            }
        )
    return sorted(rows, key=lambda row: (-row["p75_seconds"], row["test"]))


def write_csv(rows, output):
    fields = ["test", "observations", "median_seconds", "p75_seconds", "max_seconds"]
    with Path(output).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_markdown(rows, output):
    total_seconds = sum(row["median_seconds"] for row in rows)
    lines = [
        "## nf-test runtime summary",
        "",
        f"Collected **{sum(row['observations'] for row in rows)}** observations for **{len(rows)}** tests ",
        f"with **{total_seconds / 3600:.1f} runner-hours** of median execution time.",
        "",
        "| Test | Observations | Median | P75 | Maximum |",
        "|---|---:|---:|---:|---:|",
    ]
    for row in rows[:20]:
        test = row["test"].replace("|", "\\|")
        lines.append(
            f"| `{test}` | {row['observations']} | {row['median_seconds'] / 60:.1f} min | "
            f"{row['p75_seconds'] / 60:.1f} min | {row['max_seconds'] / 60:.1f} min |"
        )
    Path(output).write_text("\n".join(lines) + "\n")


def main():
    parser = argparse.ArgumentParser(description="Summarise nf-test CSV runtime reports")
    parser.add_argument("report_dir")
    parser.add_argument("--csv", required=True, dest="csv_output")
    parser.add_argument("--markdown", required=True)
    args = parser.parse_args()

    rows = summarise(read_reports(args.report_dir))
    if not rows:
        parser.error(f"no passing nf-test observations found in {args.report_dir}")
    write_csv(rows, args.csv_output)
    write_markdown(rows, args.markdown)


if __name__ == "__main__":
    main()
