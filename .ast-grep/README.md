# Nextflow structural checks

This directory owns the ast-grep CLI, Nextflow parser setup, configuration, rules
and rule tests. It does not replace Nextflow linting or pipeline execution tests.

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
The scan covers `main.nf`, `workflows/`, `subworkflows/local/` and `modules/local/`.
Imported nf-core code and configuration files remain outside this policy.

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
