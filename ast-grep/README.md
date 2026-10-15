# Nextflow structural checks

This directory owns the ast-grep CLI, Nextflow parser setup, configuration, rules
and rule tests. It does not replace Nextflow linting or pipeline execution tests.

## Setup

Install Node.js 24 with npm, a C compiler, curl, tar and either `sha256sum` or
`shasum`. macOS and Linux on arm64 and x86_64 are supported. From the repository root:

```bash
npm ci --prefix ast-grep
npm run setup --prefix ast-grep
npm run check --prefix ast-grep
```

The lockfile pins ast-grep. `setup-parser.sh` downloads a pinned, SHA-256-verified
[tree-sitter-nextflow](https://github.com/nextflow-io/tree-sitter-nextflow) source
archive and compiles its generated C parser. It retains the upstream licence in
`lib/`. Source downloads are temporary; libraries and npm dependencies are ignored
by Git. Set `CC` to a compiler executable to override `cc`. Rerun setup after
changing the parser pin or moving the checkout to another platform.

The existing prek hook runs the same check as CI. Run setup before using prek.
The scan covers `main.nf`, `workflows/`, `subworkflows/local/` and `modules/local/`.
Imported nf-core code and configuration files remain outside this policy.

## Search and add rules

Keep the configuration here and pass it explicitly:

```bash
ast-grep/node_modules/.bin/ast-grep run -c ast-grep/sgconfig.yml \
    -l nextflow -p 'if ($COND) { $$$BODY }' workflows subworkflows/local
ast-grep/node_modules/.bin/ast-grep scan -c ast-grep/sgconfig.yml --inspect summary workflows
```

Use normal `$NAME` and `$$$BODY` metavariables. `expandoChar` makes them compatible
with the Nextflow grammar. Add each policy to `rules/<id>.yml` with matching
`valid` and `invalid` snippets in `rule-tests/<id>-test.yml`. Detection-only tests
skip diagnostic snapshots. A rule with a fix must also test its rewritten output.

Warnings report findings without failing the scan. Promote an individual rule
with `--error=<id>` after reviewing its baseline. Suppress an intentional exception
with `// ast-grep-ignore: <id>` and a reason on the preceding line.

See ast-grep's [project configuration](https://ast-grep.github.io/guide/project/project-config)
and [rewriter guide](https://ast-grep.github.io/guide/rewrite/rewriter).

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
