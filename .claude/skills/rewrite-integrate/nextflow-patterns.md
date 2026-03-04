# Nextflow DSL2 Patterns for Tool Integration

Reference patterns for integrating rewritten tools into Nextflow pipelines.

## Module Binary Directory Structure

```
modules/local/<tool>/
├── main.nf
├── resources/
│   └── usr/
│       └── bin/
│           └── <binary>        # Linux x86_64 ELF binary
└── tests/                      # Optional nf-test
    └── main.nf.test
```

The binary in `resources/usr/bin/` is automatically added to `$PATH` when
`nextflow.enable.moduleBinaries = true` and Wave is enabled.

## Toggleable Process Pattern

```nextflow
// params.config or nextflow.config
params {
    use_rewritten_tool = true  // Default to rewritten version
}

// workflow
include { ORIGINAL_TOOL  } from '../modules/nf-core/tool/main'
include { REWRITTEN_TOOL } from '../modules/local/tool/main'

workflow ANALYSIS {
    take:
    ch_input
    ch_annotation

    main:
    ch_versions = Channel.empty()

    if (params.use_rewritten_tool) {
        REWRITTEN_TOOL(ch_input, ch_annotation)
        ch_results  = REWRITTEN_TOOL.out.results
        ch_versions = ch_versions.mix(REWRITTEN_TOOL.out.versions)
    } else {
        ORIGINAL_TOOL(ch_input, ch_annotation)
        ch_results  = ORIGINAL_TOOL.out.results
        ch_versions = ch_versions.mix(ORIGINAL_TOOL.out.versions)
    }

    emit:
    results  = ch_results
    versions = ch_versions
}
```

## Output Channel Compatibility

The rewritten process MUST emit channels with the same structure as the original.
If the original emits `tuple val(meta), path("*.txt")`, the rewritten must too.

Output file names and extensions must match the upstream tool exactly for MultiQC
detection. MultiQC finds tool outputs by filename patterns — if you change extensions
or naming conventions, MultiQC modules will silently skip your files.

If the rewritten tool produces files with different names, rename in the script:
```nextflow
script:
"""
<binary> ${bam} --gtf ${annotation} --outdir .

// Rename to match expected names if needed
mv rustqc_output.txt ${meta.id}.tool_output.txt
"""
```

**Glob patterns in output declarations**: Parentheses and special characters in
filenames are parsed as Groovy syntax, causing `LexerNoViableAltException`. When
upstream tools produce filenames with special characters (e.g., `counts(total).txt`),
use glob patterns in the `path()` declaration instead of literal names:
```nextflow
// BAD — parentheses break Groovy parsing
output:
tuple val(meta), path("*counts(total).txt"), emit: counts

// GOOD — glob avoids special character issues
output:
tuple val(meta), path("*total*"), emit: counts
```

## Wave Configuration

```nextflow
// nextflow.config
nextflow.enable.moduleBinaries = true

wave {
    enabled = true
    // Optional: push built containers to a registry for caching
    // build.repository = 'my-registry.io/wave-builds'
}
```

## Versions YAML Pattern

Every process should emit a `versions.yml` for pipeline provenance. Use your tool's
actual name and version in nf-core standard format (YAML with the task process name as
the top-level key):

```nextflow
output:
path "versions.yml", emit: versions

script:
"""
// ... tool execution ...

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    <binary>: \$(<binary> --version 2>&1 | head -1 | sed 's/.*version //')
END_VERSIONS
"""
```

## Preview and Lint

Always validate pipeline syntax before full runs. Run the preview after **every**
module change — it is instant and catches compilation errors immediately:

```bash
# Syntax check only — no execution, very fast
nextflow run . -profile test -preview

# Lint for common issues (if available in your Nextflow version)
nextflow lint .
```

Note: `nf-core modules install` must be run sequentially (one at a time). Parallel
installs cause git lock conflicts because each invocation modifies `modules.json`.

Add build/test artifacts to `.gitignore` immediately when starting integration work:
```
.nf-test/
test_output/
work/
```

## Resume Strategy

Always use `-resume` during integration testing:

```bash
# First run: everything executes
nextflow run . -profile test -resume --outdir results_v1

# After changing only the rewritten tool process:
# Only the changed process and its dependents re-run
nextflow run . -profile test -resume --outdir results_v2
```

This is especially valuable when upstream steps (alignment, indexing) are expensive.
Fusion-mounted previous results enable near-instant caching of unchanged steps.

## nf-test Pattern (Optional)

```groovy
// modules/local/<tool>/tests/main.nf.test
nextflow_process {
    name "Test REWRITTEN_TOOL"
    script "../main.nf"
    process "REWRITTEN_TOOL"

    test("Should produce expected outputs") {
        when {
            process {
                """
                input[0] = [
                    [ id:'test' ],
                    file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam']),
                    file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'])
                ]
                input[1] = file(params.test_data['sarscov2']['genome']['genome_gtf'])
                """
            }
        }
        then {
            assert process.success
            assert process.out.results.size() == 1
        }
    }
}
```

## Parameter Handling

When replacing a module, trace the full parameter path for any domain-specific logic.
Strandedness and paired-end logic often live in the module's `main.nf` or script
templates (reading `meta.strandedness`, `meta.single_end`), NOT in `ext.args`. If
the original module branches on metadata fields to set flags, the rewritten module
must replicate that logic — simply copying `ext.args` will miss it.

Boolean parameters that default to `false` need explicit `"default": false` in
nf-schema JSON schema files. Bare CLI `--flag` (without a value) passes the string
`"true"`, not a boolean — validate accordingly in your module.

## Common Nextflow/Groovy Pitfalls

### Conditional execution and `MissingPropertyException`

Do NOT gate process calls with `if (params.flag) { PROCESS(...) }` at the workflow
level — when the condition is false, `PROCESS.out` is never defined and downstream
references cause `MissingPropertyException`. Use `Channel.empty()` to feed the
process instead:

```nextflow
workflow ANALYSIS {
    take:
    ch_input

    main:
    ch_input_for_rewritten = params.use_rewritten_tool ? ch_input : Channel.empty()
    ch_input_for_original  = params.use_rewritten_tool ? Channel.empty() : ch_input

    REWRITTEN_TOOL(ch_input_for_rewritten)
    ORIGINAL_TOOL(ch_input_for_original)

    ch_results = REWRITTEN_TOOL.out.results.mix(ORIGINAL_TOOL.out.results)
}
```

### Disabling processes inside subworkflows with `ext.when`

To disable a process that is buried inside a subworkflow (where you cannot control
the channel wiring), use `ext.when` in config:

```nextflow
// nextflow.config
process {
    withName: '.*:SUBWORKFLOW_NAME:PROCESS_NAME' {
        ext.when = { !params.use_rewritten_tool }
    }
}
```

The process `main.nf` must check `task.ext.when` to honour this:
```nextflow
script:
def run = task.ext.when == null || task.ext.when
if (run) {
    // ... actual tool invocation ...
}
```

### String escaping in Groovy

Parentheses, brackets, and other special characters inside Nextflow `path()`
declarations are interpreted as Groovy syntax. See the glob-pattern workaround
in [Output Channel Compatibility](#output-channel-compatibility) above.
