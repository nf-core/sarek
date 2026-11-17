#!/usr/bin/env -S uv run --script
# /// script
# dependencies = []
# ///

import argparse
import re
import subprocess
import sys


def main() -> int:
    parser = argparse.ArgumentParser(description="Check added changelog lines against the current PR number.")
    parser.add_argument("--base", required=True, help="PR base commit or branch")
    parser.add_argument("--head", default="HEAD", help="PR head commit or branch")
    parser.add_argument("--pr-number", required=True, type=int)
    args = parser.parse_args()

    diff = subprocess.run(
        ["git", "diff", "--no-ext-diff", "--unified=0", f"{args.base}...{args.head}", "--", "CHANGELOG.md"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    added_lines = [line[1:] for line in diff.splitlines() if line.startswith("+") and not line.startswith("+++")]
    errors = []
    for line in added_lines:
        if re.search(r"#XXX\b|/pull/XXX\b", line, re.IGNORECASE):
            errors.append(f"Changelog PR check: replace the placeholder with the actual PR number: {line}")

    pr_link = f"[#{args.pr_number}](https://github.com/nf-core/sarek/pull/{args.pr_number})"
    if not any(line.startswith("- ") and pr_link in line for line in added_lines):
        errors.append(f"Changelog PR check: add a CHANGELOG.md entry containing {pr_link}")

    if errors:
        print("\n".join(errors), file=sys.stderr)
        return 1
    print(f"Changelog PR check: added entries reference PR #{args.pr_number} without placeholders.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
