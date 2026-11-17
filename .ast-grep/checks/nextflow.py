"""Shared helpers for Nextflow checks written with ast-grep's Python API."""

from __future__ import annotations

import os
import sys
import textwrap
from dataclasses import dataclass
from functools import cache
from pathlib import Path
from typing import Iterable, Iterator

from ast_grep_py import SgNode, SgRoot, register_dynamic_language

AST_GREP_DIR = Path(__file__).resolve().parent.parent
LANGUAGE = "nextflow"

# Directories that are never first-party pipeline code, matched as path segments.
IGNORED_DIRS = frozenset({"modules", "tests", "work", ".nf-test", ".ast-grep"})
IGNORED_SUBPATHS = (("subworkflows", "nf-core"),)


@cache
def register_nextflow() -> None:
    """Load the parser built by `mise run setup-parser` (same library as sgconfig.yml)."""
    extension = "dylib" if sys.platform == "darwin" else "so"
    library = AST_GREP_DIR / "lib" / f"libnextflow.{extension}"
    if not library.is_file():
        raise SystemExit(f"{library} is missing. Run `mise -C .ast-grep run setup-parser` first.")
    register_dynamic_language(
        {
            LANGUAGE: {
                "library_path": str(library),
                "extensions": ["nf"],
                "expando_char": "_",
            }
        }
    )


def parse(source: str) -> SgNode:
    register_nextflow()
    return SgRoot(source, LANGUAGE).root()


def is_ignored(path: Path) -> bool:
    """Match ignores against the path relative to the working directory, like the CLI."""
    parts = Path(os.path.relpath(path)).parts
    if IGNORED_DIRS.intersection(parts):
        return True
    return any(
        parts[i : i + len(sub)] == sub for sub in IGNORED_SUBPATHS for i in range(len(parts))
    )


def nextflow_files(paths: Iterable[Path]) -> Iterator[Path]:
    """Expand files and directories into `.nf` files, skipping ignored paths."""
    for path in paths:
        candidates = sorted(path.rglob("*.nf")) if path.is_dir() else [path]
        for candidate in candidates:
            if candidate.suffix == ".nf" and not is_ignored(candidate):
                yield candidate


def is_suppressed(node: SgNode, rule_id: str) -> bool:
    """Honour `// ast-grep-ignore` or `// ast-grep-ignore: <id>[, <id>]` on the line above."""
    previous = node.prev()
    if previous is None or previous.kind() != "line_comment":
        return False
    if previous.range().end.line != node.range().start.line - 1:
        return False
    text = previous.text().removeprefix("//").strip()
    if not text.startswith("ast-grep-ignore"):
        return False
    _, _, ids = text.partition(":")
    return not ids.strip() or rule_id in {i.strip() for i in ids.split(",")}


@dataclass(frozen=True)
class Finding:
    rule_id: str
    path: Path
    node: SgNode
    message: str

    def format(self) -> str:
        start = self.node.range().start
        source = textwrap.dedent(" " * start.column + self.node.text())
        snippet = "\n".join(f"    | {line}" for line in source.splitlines())
        return (
            f"warning[{self.rule_id}]: {self.message}\n"
            f"  --> {self.path}:{start.line + 1}:{start.column + 1}\n{snippet}"
        )
