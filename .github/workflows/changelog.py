#! /usr/bin/env python3
"""
Taken from https://github.com/MultiQC/MultiQC/blob/main/.github/workflows/changelog.py and updated for nf-core

To be called by a CI action. Assumes the following environment variables are set:
PR_TITLE, PR_NUMBER, GITHUB_WORKSPACE.

Adds a line into the CHANGELOG.md:
* Looks for the section to add the line to, based on the PR title, e.g. `Added:`, `Changed:`.
* All other change will go under the "## dev" section.
* If an entry for the PR is already added, it will not run.

Other assumptions:
- CHANGELOG.md has a running section for an ongoing "dev" version
(i.e. titled "## dev").
"""

import os
import re
import sys
from pathlib import Path
from typing import List, Tuple

REPO_URL = "https://github.com/nf-core/tools"

# Assumes the environment is set by the GitHub action.
pr_title = os.environ["PR_TITLE"]
pr_number = os.environ["PR_NUMBER"]
comment = os.environ.get("COMMENT", "")
workspace_path = Path(os.environ.get("GITHUB_WORKSPACE", ""))

assert pr_title, pr_title
assert pr_number, pr_number

# Trim the PR number added when GitHub squashes commits, e.g. "Added: Updated (#2026)"
pr_title = pr_title.removesuffix(f" (#{pr_number})")  # type: ignore

changelog_path = workspace_path / "CHANGELOG.md"

if any(
    line in pr_title.lower()
    for line in [
        "skip changelog",
        "skip change log",
        "no changelog",
        "no change log",
        "bump version",
    ]
):
    print("Skipping changelog update")
    sys.exit(0)


def _determine_change_type(pr_title) -> Tuple[str, str]:
    """
    Determine the type of the PR: Added, Changed, Fixed, or Removed
    Returns a tuple of the section name and the module info.
    """
    sections = {
        "Added": "### Added",
        "Changed": "### Changed",
        "Fixed": "### Fixed",
        "Removed": "### Removed",
    }
    current_section_header = "## dev"
    current_section = "dev"

    # Check if the PR in any of the sections.
    for section, section_header in sections.items():
        # check if the PR title contains any of the section headers, with some loose matching, e.g. removing plural and suffixes
        if re.sub(r"s$", "", section.lower().replace("ing", "")) in pr_title.lower():
            current_section_header = section_header
            current_section = section
    print(f"Detected section: {current_section}")
    return current_section, current_section_header


# Determine the type of the PR
section, section_header = _determine_change_type(pr_title)

# Remove section indicator from the PR title.
pr_title = re.sub(rf"{section}:[\s]*", "", pr_title, flags=re.IGNORECASE)

# Prepare the change log entry.
pr_link = f"([#{pr_number}]({REPO_URL}/pull/{pr_number}))"

# Handle manual changelog entries through comments.
if comment := comment.removeprefix("@nf-core-bot changelog").strip():  # type: ignore
    print(f"Adding manual changelog entry: {comment}")
    pr_title = comment
new_lines = [
    f"- {pr_link} - {pr_title}\n",
]
print(f"Adding new lines into section '{section}':\n" + "".join(new_lines))

# Finally, updating the changelog.
# Read the current changelog lines. We will print them back as is, except for one new
# entry, corresponding to this new PR.
with changelog_path.open("r") as f:
    orig_lines = f.readlines()
updated_lines: List[str] = []


def _skip_existing_entry_for_this_pr(line: str, same_section: bool = True) -> str:
    if line.strip().endswith(pr_link):
        print(f"Found existing entry for this pull request #{pr_number}:")
        existing_lines = [line]
        if new_lines and new_lines == existing_lines and same_section:
            print(
                f"Found existing identical entry for this pull request #{pr_number} in the same section:"
            )
            print("".join(existing_lines))
            sys.exit(0)  # Just leaving the CHANGELOG intact
        else:
            print(
                f"Found existing entry for this pull request #{pr_number}. It will be replaced and/or moved to proper section"
            )
            print("".join(existing_lines))
            for _ in range(len(existing_lines)):
                try:
                    line = orig_lines.pop(0)
                except IndexError:
                    break
    return line


# Find the next line in the change log that matches the pattern "## dev"
# If it doesn't exist, exist with code 1 (let's assume that a new section is added
# manually or by CI when a release is pushed).
# Else, find the next line that matches the `section` variable, and insert a new line
# under it (we also assume that section headers are added already).
inside_version_dev = False
already_added_entry = False
while orig_lines:
    line = orig_lines.pop(0)

    # If the line already contains a link to the PR, don't add it again.
    line = _skip_existing_entry_for_this_pr(line, same_section=False)

    if (
        line.startswith("## ") and not line.strip() == "# nf-core/sarek: Changelog"
    ):  # Version header, e.g. "## 2.12dev"
        print(f"Found version header: {line.strip()}")
        updated_lines.append(line)

        # Parse version from the line `## 3.4.4` or
        # `## [3.4.4](https://github.com/nf-core/sarek/releases/tag/3.4.4) - Ruopsokjåkhå` ...
        if not (m := re.match(r".*(\d+\.\d+.\d*(dev)?).*", line)):
            print(f"Cannot parse version from line {line.strip()}.", file=sys.stderr)
            sys.exit(1)
        version = m.group(1)
        print(f"Found version: {version}")

        if not inside_version_dev:
            if not version.endswith("dev"):
                print(
                    "Can't find a 'dev' version section in the changelog. Make sure "
                    "it's created, and all the required sections, e.g. `### Template` are created under it .",
                    file=sys.stderr,
                )
                sys.exit(1)
            inside_version_dev = True
        else:
            if version.endswith("dev"):
                print(
                    f"Found another 'dev' version section in the changelog, make"
                    f"sure to change it to a 'release' stable version tag. "
                    f"Line: {line.strip()}",
                    file=sys.stderr,
                )
                sys.exit(1)
            # We are past the dev version, so just add back the rest of the lines and break.
            while orig_lines:
                line = orig_lines.pop(0)
                line = _skip_existing_entry_for_this_pr(line, same_section=False)
                if line:
                    updated_lines.append(line)
            break
        continue
    print(f"Found line: {line.strip()}")
    print(f"inside_version_dev: {inside_version_dev}")
    print(f"section_header: {section_header}")
    if inside_version_dev and line.lower().startswith(
        section_header.lower()
    ):  # Section of interest header
        print(f"Found section header: {line.strip()}")
        if already_added_entry:
            print(
                f"Already added new lines into section {section}, is the section duplicated?",
                file=sys.stderr,
            )
            sys.exit(1)
        updated_lines.append(line)
        # Collecting lines until the next section.
        section_lines: List[str] = []
        while True:
            line = orig_lines.pop(0)
            if line.startswith("#"):
                print(f"Found the next section header: {line.strip()}")
                # Found the next section header, so need to put all the lines we collected.
                updated_lines.append("\n")
                _updated_lines = [_l for _l in section_lines + new_lines if _l.strip()]
                updated_lines.extend(_updated_lines)
                updated_lines.append("\n")
                if new_lines:
                    print(
                        f"Updated {changelog_path} section '{section}' with lines:\n"
                        + "".join(new_lines)
                    )
                else:
                    print(
                        f"Removed existing entry from {changelog_path} section '{section}'"
                    )
                already_added_entry = True
                # Pushing back the next section header line
                orig_lines.insert(0, line)
                break
            # If the line already contains a link to the PR, don't add it again.
            line = _skip_existing_entry_for_this_pr(line, same_section=True)
            section_lines.append(line)

    else:
        updated_lines.append(line)


def collapse_newlines(lines: List[str]) -> List[str]:
    updated = []
    for idx in range(len(lines)):
        if idx != 0 and not lines[idx].strip() and not lines[idx - 1].strip():
            continue
        updated.append(lines[idx])
    return updated


updated_lines = collapse_newlines(updated_lines)

# Finally, writing the updated lines back.
with changelog_path.open("w") as f:
    f.writelines(updated_lines)
