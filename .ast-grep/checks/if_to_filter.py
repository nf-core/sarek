"""Warn when a workflow `if` guard could become a filter on the driving input channel.

Matches a braced `if` without `else` that sits directly in a workflow body or
`main:` section and contains exactly one statement (comments aside): an
unqualified call to an uppercase name, for example `ALIGN(reads, reference)`.

Uppercase naming suggests a process or subworkflow but does not prove it. The
check only reports; converting a guard to a filter needs human review.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Iterable, Iterator, Optional

from ast_grep_py import SgNode

from checks.nextflow import Finding, is_suppressed, nextflow_files, parse

RULE_ID = "if-to-filter"
PROCESS_NAME = re.compile(r"^[A-Z][A-Z0-9_]*$")
COMMENTS = frozenset({"line_comment", "block_comment"})
WORKFLOW_SCOPES = frozenset({"workflow_definition", "main_section"})
NOTE = """\
Review whether the condition can move into a filter on the execution-driving
queue channel. Preserve other inputs and review output wiring and empty-input
behaviour. This is a naming-based suggestion, not a resolved process call.
No automatic fix: channel types and subworkflow side effects require review."""


def guarded_call(guard: SgNode) -> Optional[SgNode]:
    """Return the single uppercase call gated by `guard`, or None if it does not qualify."""
    if guard.field("alternative") is not None:
        return None
    parent = guard.parent()
    if parent is None or parent.kind() not in WORKFLOW_SCOPES:
        return None
    block = guard.field("consequence")
    if block is None or block.kind() != "block":
        return None

    statements = [c for c in block.children() if c.is_named() and c.kind() not in COMMENTS]
    if len(statements) != 1 or statements[0].kind() != "expression_statement":
        return None
    call = statements[0].child(0)
    if call is None or call.kind() != "call_expression":
        return None
    name = call.child(0)
    if name is None or name.kind() != "identifier" or not PROCESS_NAME.match(name.text()):
        return None
    return call


def check_source(source: str, path: Path = Path("<stdin>")) -> Iterator[Finding]:
    for guard in parse(source).find_all(kind="if_statement"):
        call = guarded_call(guard)
        if call is None or is_suppressed(guard, RULE_ID):
            continue
        process = call.child(0).text()
        message = f"Prefer filtering the driving input channel to conditionally invoking {process}."
        yield Finding(RULE_ID, path, guard, message)


def check_paths(paths: Iterable[Path]) -> list[Finding]:
    return [
        finding
        for path in nextflow_files(paths)
        for finding in check_source(path.read_text(), path)
    ]


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("paths", nargs="+", type=Path, help="Nextflow files or directories")
    parser.add_argument(
        "--error", action="store_true", help="exit non-zero when any guard is reported"
    )
    args = parser.parse_args(argv)

    missing = [str(p) for p in args.paths if not p.exists()]
    if missing:
        parser.error(f"no such path: {', '.join(missing)}")

    findings = check_paths(args.paths)
    for finding in findings:
        print(finding.format(), end="\n\n")
    if findings:
        print(f"note[{RULE_ID}]:\n{NOTE}\n")
    print(f"{RULE_ID}: {len(findings)} warning(s)")
    return 1 if args.error and findings else 0


if __name__ == "__main__":
    sys.exit(main())
