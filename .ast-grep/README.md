# Structural checks

This directory owns the ast-grep CLI, Nextflow parser setup, configuration, rules
and rule tests for Nextflow and the Markdown changelog. It does not replace
Nextflow linting or pipeline execution tests.

## Setup

Install [mise](https://mise.jdx.dev/walkthrough.html), a C compiler, curl, tar and either `sha256sum` or
`shasum`. macOS and Linux on arm64 and x86_64 are supported. From the repository root:

```bash
mise trust .ast-grep/mise.toml
mise -C .ast-grep run setup-parser
mise -C .ast-grep run check
```

Mise installs the pinned native ast-grep binary. Node.js and npm are not required.
The file task `mise-tasks/setup-parser` downloads a pinned, SHA-256-verified
[tree-sitter-nextflow](https://github.com/nextflow-io/tree-sitter-nextflow) source
archive and compiles its generated C parser. It retains the upstream licence in
`lib/`. Source downloads are temporary; libraries are ignored by Git.
The task also smoke-tests the parser. `setup` is an alias for `setup-parser`.
Set `CC` to a compiler executable to override `cc`. Rerun setup-parser after
changing the parser pin or moving the checkout to another platform.

The existing prek hook runs the same check as CI. Run setup before using prek.
The scan covers `main.nf`, `workflows/`, `subworkflows/local/`, `modules/local/` and
`CHANGELOG.md`. Imported nf-core code and configuration files remain outside this
policy. Markdown support is built into the pinned ast-grep binary.

## Search and add rules

Keep the configuration here and pass it explicitly:

```bash
mise -C .ast-grep exec -- ast-grep run -c sgconfig.yml \
    -l nextflow -p 'if ($COND) { $$$BODY }' ../workflows ../subworkflows/local
mise -C .ast-grep exec -- ast-grep scan -c sgconfig.yml --inspect summary ../workflows
```

Use normal `$NAME` and `$$$BODY` metavariables. `expandoChar` makes them compatible
with the Nextflow grammar. Add each policy to `rules/<id>.yml` with matching
`valid` and `invalid` snippets in `rule-tests/<id>-test.yml`. Generate snapshots with
`mise -C .ast-grep exec -- ast-grep test -c sgconfig.yml --update-all`, then review
and commit them. The `check` task verifies both matching and snapshot results.
A rule with a fix must also test its rewritten output.

Warnings report findings without failing the scan. Promote an individual rule
with `--error=<id>` after reviewing its baseline. Suppress an intentional exception
with `// ast-grep-ignore: <id>` and a reason on the preceding line.

See ast-grep's [project configuration](https://ast-grep.github.io/guide/project/project-config)
and [rewriter guide](https://ast-grep.github.io/guide/rewrite/rewriter).

## Changelog PR references

`changelog-pr-placeholder` reports placeholder PR links inside Markdown list items
in `CHANGELOG.md`. Fenced code examples and other Markdown files are excluded.
The local hook runs its native YAML tests and scans the changelog.

CI also runs `check-changelog-pr`. This task extracts added changelog lines from
the PR's merge-base diff, then uses ast-grep Markdown rules to reject placeholders
and require a list item linking to the current PR number. Both the link label and
URL must use that number. The shell task supplies Git and PR metadata; matching
and diagnostics come from ast-grep. Unchanged historical entries are excluded
from this PR-diff check.

```bash
BASE_SHA=origin/dev PR_NUMBER=2424 mise -C .ast-grep run check-changelog-pr
```

Set `BASE_SHA` to the actual parent branch or commit for a stacked PR. Set
`HEAD_SHA` to inspect a committed revision other than `HEAD`. The normal `check`
task also runs temporary-Git-repository tests of the PR-diff task.

Open a draft PR to obtain its number before committing a new changelog entry.
Push the numbered entry before requesting review.

## Convert an if guard to a filter

`if-to-filter` warns about a braced `if` containing a single uppercase call directly
in a workflow's implicit body or `main:` section. Comments are allowed. It excludes
`else`, nested guards, multi-statement bodies, assignments, closures, functions and
qualified calls. Uppercase naming suggests a process or subworkflow but does not
resolve its definition. Modules, imported subworkflows and test/work directories
are excluded even when scanning a larger tree.

For a process driven by a queue channel, review this manual conversion:

```nextflow
// Before
if (params.align) {
    ALIGN(reads, reference)
}

// After, when reads is the execution-driving queue channel
ALIGN(reads.filter { _read -> params.align }, reference)
```

There is no automatic fix. A filter evaluates its predicate per item, whereas the
guard evaluates it once. Confirm that the condition is stable and side-effect-free,
and that the gated input is a queue channel. Review other inputs, output wiring and
empty-input behaviour. Invoking a subworkflow with an empty input can still execute
independent work inside it. Test enabled, disabled and empty-input cases before
applying the conversion. Multi-call and `else` blocks need a separate design review.
