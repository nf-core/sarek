#!/usr/bin/env python3

import argparse
import csv
import json
import re
from pathlib import Path


ANSI_ESCAPE = re.compile(r"\x1b\[[0-9;]*m")
HASH_PATTERN = re.compile(r"\bTest\b.*\[([0-9a-f]{8})\]")


def normalise_path(filename):
    path = Path(filename)
    parts = path.parts
    if "tests" in parts:
        return Path(*parts[parts.index("tests") :]).as_posix()
    return path.as_posix()


def read_discovery(discovery_csv, dry_run_log):
    with Path(discovery_csv).open(newline="") as handle:
        rows = list(csv.DictReader(handle))

    hashes = []
    for line in Path(dry_run_log).read_text().splitlines():
        match = HASH_PATTERN.search(ANSI_ESCAPE.sub("", line))
        if match:
            hashes.append(match.group(1))

    if len(rows) != len(hashes):
        raise ValueError(f"discovery CSV contains {len(rows)} tests but dry-run log contains {len(hashes)} hashes")

    tests = []
    for row, test_hash in zip(rows, hashes, strict=True):
        path = normalise_path(row["filename"])
        tests.append(
            {
                "key": f"{path}::{row['test']}",
                "selector": f"{path}@{test_hash}",
            }
        )
    return tests


def read_weights(weights_csv):
    with Path(weights_csv).open(newline="") as handle:
        return {row["test"]: float(row["p75_seconds"]) for row in csv.DictReader(handle)}


def assign_shards(tests, weights, shard_count, default_weight=None):
    known_weights = sorted(weights.values())
    if default_weight is None:
        default_weight = known_weights[len(known_weights) * 3 // 4] if known_weights else 120.0

    weighted_tests = [
        (weights.get(test["key"], default_weight), test["selector"])
        for test in tests
    ]
    weighted_tests.sort(key=lambda item: (-item[0], item[1]))

    shards = [{"index": index + 1, "estimated_seconds": 0.0, "selectors": []} for index in range(shard_count)]
    for weight, selector in weighted_tests:
        shard = min(shards, key=lambda item: (item["estimated_seconds"], item["index"]))
        shard["selectors"].append(selector)
        shard["estimated_seconds"] += weight

    for shard in shards:
        shard["estimated_seconds"] = round(shard["estimated_seconds"], 3)
    return shards


def validate_assignments(tests, shards):
    expected = sorted(test["selector"] for test in tests)
    assigned = sorted(selector for shard in shards for selector in shard["selectors"])
    if assigned != expected:
        raise ValueError("weighted shard assignments do not match discovered tests")


def main():
    parser = argparse.ArgumentParser(description="Create duration-weighted nf-test shards")
    parser.add_argument("--discovery-csv", required=True)
    parser.add_argument("--dry-run-log", required=True)
    parser.add_argument("--weights", required=True)
    parser.add_argument("--shards", required=True, type=int)
    args = parser.parse_args()

    tests = read_discovery(args.discovery_csv, args.dry_run_log)
    shards = assign_shards(tests, read_weights(args.weights), args.shards)
    validate_assignments(tests, shards)
    print(json.dumps(shards, separators=(",", ":")))


if __name__ == "__main__":
    main()
