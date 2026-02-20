# nf-core/sarek Developer Guidelines

This document provides comprehensive guidelines for contributing to the nf-core/sarek pipeline. These guidelines are designed for both human developers and AI agents.

## Table of Contents

- [Contributing Principles](#contributing-principles)
- [Git Workflow](#git-workflow)
- [Codebase Architecture](#codebase-architecture)
- [Code Style](#code-style)
- [Channel Operations and Gotchas](#channel-operations-and-gotchas)
- [Meta Map Handling](#meta-map-handling)
- [Modules](#modules)
- [Subworkflows](#subworkflows)
- [Configuration](#configuration)
- [Testing](#testing)
- [Documentation](#documentation)
- [Metro Map Updates](#metro-map-updates)
- [PR Checklist](#pr-checklist)

---

## Contributing Principles

- **Read files before editing** — understand existing code before making changes
- Keep fixes **minimal and focused** — don't refactor surrounding code
- Don't add docstrings, comments, or type annotations to unchanged code
- Don't add error handling or validation beyond what's needed
- Don't over-engineer: no premature abstractions, no feature flags
- When unsure about scope or approach, ask rather than guess

---

## Git Workflow

- **Always branch off `origin/dev`**, never master
- Branch naming: `fix/issue-XXXX` or `feat/issue-XXXX`
- PRs target the `dev` branch
- Never force push, never amend published commits without asking
- Commit messages should be descriptive and include the issue reference

---

## Codebase Architecture

Sarek follows a hierarchical, modular architecture:

```
Modules (atomic processes) → Subworkflows (composed modules) → Workflow (orchestration)
```

**Key design principles:**

- Separation of concerns between processing steps
- Reusable components through nf-core modules ecosystem
- Configuration-driven behavior via `ext.*` directives
- Comprehensive testing with nf-test

### Directory Structure

```
sarek/
├── main.nf                     # Pipeline entry point
├── nextflow.config             # Main configuration
├── nextflow_schema.json        # Parameter schema (JSON Schema)
├── modules.json                # nf-core module tracking
├── modules/
│   ├── local/                  # Pipeline-specific modules
│   └── nf-core/                # Imported nf-core modules (120+)
├── subworkflows/
│   ├── local/                  # Pipeline-specific subworkflows (66)
│   └── nf-core/                # Imported nf-core subworkflows (6)
├── workflows/sarek/main.nf     # Main workflow orchestration
├── conf/
│   ├── base.config             # Default resource allocations
│   ├── modules/                # Module-specific configurations (41 files)
│   └── test/                   # Test configurations (21 files)
├── tests/                      # nf-test test files (57 tests)
├── docs/                       # Documentation
└── assets/                     # MultiQC config, samplesheets, etc.
```

---

## Code Style

### Harshil Alignment

Use "Harshil alignment" for include statements - align the closing braces to improve readability:

```groovy
// CORRECT - Harshil alignment
include { paramsSummaryMap                                  } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                              } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                            } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                            } from '../../subworkflows/local/utils_nfcore_sarek_pipeline'

// CORRECT - With aliases
include { BAM_CONVERT_SAMTOOLS as CONVERT_FASTQ_INPUT       } from '../../subworkflows/local/bam_convert_samtools'
include { SPRING_DECOMPRESS as SPRING_DECOMPRESS_TO_R1_FQ   } from '../../modules/nf-core/spring/decompress'
include { SPRING_DECOMPRESS as SPRING_DECOMPRESS_TO_R2_FQ   } from '../../modules/nf-core/spring/decompress'

// INCORRECT - No alignment
include { paramsSummaryMap } from 'plugin/nf-schema'
include { paramsSummaryMultiqc } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
```

### Harshil Alignment in Take/Emit Blocks

Also apply alignment to `take:` and `emit:` blocks:

```groovy
take:
cram          // channel: [mandatory] [ meta, cram, crai ]
dict          // channel: [optional]  [ meta, dict ]
fasta         // channel: [mandatory] [ fasta ]
fasta_fai     // channel: [mandatory] [ fasta_fai ]
intervals     // channel: [mandatory] [ interval.bed.gz, interval.bed.gz.tbi, num_intervals ]

emit:
vcf_ann       // channel: [ val(meta), vcf.gz, vcf.gz.tbi ]
tab_ann
json_ann
reports       //    path: *.html
versions      //    path: versions.yml
```

### Channel Naming Conventions

```groovy
// Initial process output channel
ch_output_from_<process>

// Intermediate/terminal channels
ch_<previousprocess>_for_<nextprocess>

// Example
ch_bam_from_markduplicates
ch_markduplicates_for_baserecalibrator
```

### Topic Channels

We are migrating to **Nextflow topic channels** where possible. Topics allow processes and subworkflows to publish to a named topic without explicit channel wiring.

**When touching code in a PR, migrate `versions.mix(...)` to topic channels:**

```groovy
// OLD - Explicit version collection
versions = versions.mix(TOOL_A.out.versions)
versions = versions.mix(TOOL_B.out.versions)

// NEW - Use topic channel
// In the process/module, publish versions to a topic:
//   output:
//   path "versions.yml", topic: 'versions'
//
// In the workflow, collect from the topic:
//   ch_versions = Channel.topic('versions')
```

**Where to use topics:**

- Version collection (`versions` topic) — primary migration target
- Any cross-subworkflow channel passing that doesn't depend on meta-keyed joins

**Where NOT to use topics:**

- Channels that need `join`, `groupTuple`, or other keyed operations
- Channels where ordering or pairing matters

### General Style

- Use 4-space indentation
- Put channel operations on separate lines for readability
- Add comments for complex logic
- Use descriptive variable names

### Strict Syntax Mode

**When touching any code in a PR, you must update it to use strict Nextflow syntax.** This ensures gradual modernization of the codebase.

#### Required Changes When Modifying Code

1. **Use explicit `it` variable or named parameters in closures:**

   ```groovy
   // CORRECT - Explicit named parameters
   .map { meta, vcf -> [meta, vcf] }

   // CORRECT - Explicit `it` when single parameter
   .map { it -> it.baseName }

   // DEPRECATED - Implicit `it`
   .map { it.baseName }
   ```

2. **Use parentheses for all method calls:**

   ```groovy
   // CORRECT
   channel.map({ meta, file -> [meta, file] })

   // ALSO CORRECT (trailing closure)
   channel.map { meta, file -> [meta, file] }
   ```

3. **Explicit type declarations where applicable:**

   ```groovy
   // CORRECT
   String prefix = "${meta.id}"
   List<String> args = []

   // AVOID in new code
   def prefix = "${meta.id}"
   ```

4. **Use underscore prefix for unused/dropped variables:**

   The underscore prefix convention clearly indicates which variables from a closure are intentionally not used in the output. This makes code review easier and prevents confusion about whether a variable was accidentally omitted.

   ```groovy
   // CORRECT - Underscore prefix shows vcf is intentionally dropped
   .map { meta, _vcf, tbi -> [meta, tbi] }

   // CORRECT - Multiple dropped variables
   .map { meta, _vcf, _tbi, file -> [meta, file] }

   // CORRECT - In join operations
   .join(other_channel, failOnDuplicate: true, failOnMismatch: true)
       .map { meta, file1, _file2 -> [meta, file1] }

   // CORRECT - When extracting from complex structures
   VCF_ANNOTATE_SNPEFF.out.vcf_tbi.map { meta, vcf_, _tbi -> [meta, vcf_, []] }

   // INCORRECT - Unclear which variables are intentionally unused
   .map { meta, vcf, tbi -> [meta, tbi] }
   ```

   **When to use underscore prefix:**
   - Variable is received but not included in output
   - Variable is needed for destructuring but value is discarded
   - Makes intent clear during code review

---

## Channel Operations and Gotchas

### Join Operations - ALWAYS Use `failOnDuplicate` and `failOnMismatch`

When joining channels, ALWAYS specify `failOnDuplicate: true, failOnMismatch: true` to catch bugs early:

```groovy
// CORRECT - Will fail fast if there are issues
vcf_tbi = vcf.join(tbi, failOnDuplicate: true, failOnMismatch: true)

// INCORRECT - Silent failures can cause subtle bugs
vcf_tbi = vcf.join(tbi)
```

Use `remainder: true` only when intentionally handling unmatched items:

```groovy
// When some items may not have matches (intentional)
all_unmapped_bam = SAMTOOLS_VIEW_UNMAP_UNMAP.out.bam
    .join(SAMTOOLS_VIEW_UNMAP_MAP.out.bam, failOnDuplicate: true, remainder: true)
    .join(SAMTOOLS_VIEW_MAP_UNMAP.out.bam, failOnDuplicate: true, remainder: true)
```

### Branch Operations

Use `branch` to split channels based on conditions:

```groovy
vcf_out = STRELKA_SINGLE.out.vcf.branch{
    // Use meta.num_intervals to assess number of intervals
    intervals:    it[0].num_intervals > 1
    no_intervals: it[0].num_intervals <= 1
}

// Access branches
vcf_out.intervals      // Items where num_intervals > 1
vcf_out.no_intervals   // Items where num_intervals <= 1
```

### GroupTuple - Use `groupKey` for Performance

When using `groupTuple`, use `groupKey` with known size to avoid blocking:

```groovy
// CORRECT - Non-blocking when size is known
vcf_to_merge = vcf_out.intervals
    .map{ meta, vcf -> [ groupKey(meta, meta.num_intervals), vcf ]}
    .groupTuple()

// NOTE: Without groupKey and size, groupTuple is a blocking operation
// This can cause pipeline hangs if the expected number of items varies
```

### Strelka Special Case - SNV and Indel VCFs

Strelka produces TWO VCF files (SNVs and Indels) that need special handling:

```groovy
// Strelka somatic outputs need to be concatenated before consensus calling
ch_vcfs = vcfs.branch{ meta, vcf, tbi ->
    strelka_somatic: meta.variantcaller == 'strelka' && meta.status == '1'
    other: true
}

// Concatenate the two strelka VCFs (SNPs and indels) using groupTuple(size: 2)
BCFTOOLS_CONCAT(ch_vcfs.strelka_somatic.groupTuple(size: 2))
```

### Combine vs Join

- Use `join` when combining channels by a key (meta map)
- Use `combine` when creating cartesian product (e.g., sample x intervals)

```groovy
// Join by meta key
vcf_tbi = vcf.join(tbi, failOnDuplicate: true, failOnMismatch: true)

// Combine all samples with all intervals (cartesian product)
cram_intervals = cram.combine(intervals)
```

### Controlling Flow with Channel Operations (Preferred)

**Nextflow is a dataflow language.** Prefer channel operations over `if` statements to control which processes run:

```groovy
// BEST - Use filter to control what enters a process
input_channel
    .filter { meta, _file -> params.tools?.split(',')?.contains('toolname') }
    .set { ch_for_tool }

TOOL_PROCESS(ch_for_tool)

// BEST - Use branch for multiple conditional paths
input_channel.branch { meta, file ->
    tool_a: params.tools?.split(',')?.contains('tool_a')
    tool_b: params.tools?.split(',')?.contains('tool_b')
    other: true
}.set { ch_branched }

TOOL_A(ch_branched.tool_a)
TOOL_B(ch_branched.tool_b)

// AVOID - if statements for flow control (use only when channel ops aren't suitable)
if (params.run_tool) {
    TOOL_PROCESS(input_channel)
}
```

**Benefits of channel operations:**

- More idiomatic Nextflow - data drives execution
- Better composability and testability
- Clearer dataflow visualization
- Avoids caching issues when conditions change

---

## Meta Map Handling

### Adding Fields to Meta

Use `meta + [key: value]` syntax:

```groovy
// Add single field
meta = meta + [id: meta.sample]

// Add multiple fields
meta = meta + [id: "${meta.sample}-${meta.lane}".toString(), data_type: "fastq_gz", num_lanes: num_lanes.toInteger()]

// In map operation
.map{ meta, vcf -> [ meta + [ variantcaller:'strelka' ], vcf ] }
```

### Removing Fields from Meta - Use `subMap`

Use `meta - meta.subMap('field')` to remove fields:

```groovy
// Remove single field
.map{ meta, vcf -> [ meta - meta.subMap('num_intervals'), vcf ] }

// Remove multiple fields
.map{ meta, vcf, tbi ->
    [meta - meta.subMap('variantcaller', 'contamination', 'filename'), vcf, tbi]
}

// Add and remove in one operation
.map{ meta, vcf -> [ meta - meta.subMap('num_intervals') + [ variantcaller:'strelka' ], vcf ] }
```

### Accessing Meta Fields

```groovy
// In map closures
.map{ meta, file -> [meta.sample, file] }

// In branch conditions
.branch{ meta, vcf ->
    intervals:    meta.num_intervals > 1
    no_intervals: meta.num_intervals <= 1
}

// Getting subset of meta
[meta.patient, meta.subMap('sample', 'status')]
```

### Common Meta Fields in Sarek

| Field                | Description                                   |
| -------------------- | --------------------------------------------- |
| `meta.patient`       | Patient identifier                            |
| `meta.sample`        | Sample identifier                             |
| `meta.status`        | 0 = normal, 1 = tumor                         |
| `meta.lane`          | Sequencing lane                               |
| `meta.id`            | Unique identifier (often `${sample}-${lane}`) |
| `meta.data_type`     | Input type: `fastq_gz`, `bam`, `cram`         |
| `meta.num_intervals` | Number of intervals for scatter/gather        |
| `meta.variantcaller` | Name of variant caller                        |
| `meta.num_lanes`     | Total number of lanes for sample              |

---

## Modules

### DEPRECATED: The `ext.when` Clause Pattern

> **DEPRECATED:** The `ext.when` clause pattern is deprecated and should NOT be used in new code. Existing code using this pattern should be refactored when touched in a PR.

You may see comments in older subworkflow files like:

```groovy
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run
```

**Do not follow this pattern for new code.** Instead, prefer channel operations to control dataflow:

```groovy
// BEST - Use channel operations (filter, branch) to control dataflow
input_channel
    .filter { meta, _file -> params.tools?.split(',')?.contains('toolname') }
    .set { ch_for_tool }

TOOL_PROCESS(ch_for_tool)

// ACCEPTABLE - When channel operations aren't suitable, use explicit conditional
if (tools && tools.split(',').contains('toolname')) {
    TOOL_PROCESS(input_channel)
    versions = versions.mix(TOOL_PROCESS.out.versions)
}

// DEPRECATED - Using ext.when in config
// withName: 'TOOL_PROCESS' {
//     ext.when = { params.tools && params.tools.split(',').contains('toolname') }
// }
```

**Why channel operations are preferred:**

- Nextflow is a dataflow language - let the data drive execution
- Channel operations are more composable and testable
- Avoids issues with process caching when conditions change
- Makes the pipeline logic more explicit and traceable

### Remapping Channels for Module Input

When a module expects different input structure, remap in the call:

```groovy
// Remap channel to match module/subworkflow input signature
BAM_VARIANT_CALLING_CNVKIT(
    cram.map{ meta, cram, crai -> [ meta, [], cram ] },
    fasta,
    fasta_fai,
    intervals_bed_combined.map{it -> it ? [[id:it[0].baseName], it]: [[id:'no_intervals'], []]},
    params.cnvkit_reference ? cnvkit_reference.map{ it -> [[id:it[0].baseName], it] } : [[:],[]]
)
```

### Module Memory Requirements

Some modules have specific memory requirements noted in comments:

```groovy
// In modules/nf-core/bwa/index/main.nf:
// NOTE requires 5.37N memory where N is the size of the database

// In modules/nf-core/bwamem2/index/main.nf:
// NOTE Requires 28N GB memory where N is the size of the reference sequence, floor of 280M
```

### Adding/Updating nf-core Modules

```bash
# Install a new module
nf-core modules install <tool>/<subcommand>

# Update an existing module
nf-core modules update <tool>/<subcommand>

# List installed modules
nf-core modules list local
```

---

## Subworkflows

### Subworkflow Naming Patterns

| Category        | Naming Pattern          | Examples                                  |
| --------------- | ----------------------- | ----------------------------------------- |
| Alignment       | `fq_align_*`            | `fq_align_bwamem`, `fq_align_bwamem2`     |
| BAM processing  | `bam_*`                 | `bam_markduplicates`, `bam_applybqsr`     |
| Variant calling | `bam_variant_calling_*` | `bam_variant_calling_germline_all`        |
| VCF processing  | `vcf_*`                 | `vcf_annotate_all`, `vcf_concat_variants` |
| Preparation     | `prepare_*`             | `prepare_genome`, `prepare_intervals`     |

### Subworkflow Structure

```groovy
//
// DESCRIPTION OF SUBWORKFLOW
//

include { MODULE_A                    } from '../../../modules/nf-core/module_a'
include { MODULE_B                    } from '../../../modules/nf-core/module_b'
include { MODULE_B as MODULE_B_ALIAS  } from '../../../modules/nf-core/module_b'

workflow SUBWORKFLOW_NAME {
    take:
    input_channel    // channel: [mandatory] [ meta, file ]
    other_inputs     // channel: [optional]  description

    main:
    versions = Channel.empty()

    // Initialize output channels
    output_a = Channel.empty()
    output_b = Channel.empty()

    // PREFERRED: Use channel operations to control dataflow
    ch_for_module_a = input_channel
        .filter { meta, _file -> meta.run_module_a }

    MODULE_A(ch_for_module_a)
    versions = versions.mix(MODULE_A.out.versions)

    MODULE_B(MODULE_A.out.result)
    versions = versions.mix(MODULE_B.out.versions)

    emit:
    result   = MODULE_B.out.result  // channel: [ val(meta), file ]
    versions                        // channel: versions.yml
}
```

### Scatter-Gather Pattern

Common pattern for parallelizing over intervals:

```groovy
// Combine samples with intervals for scatter strategy
cram_intervals = cram.combine(intervals)
    // Move num_intervals to meta map for later grouping
    .map{ meta, cram, crai, intervals, intervals_index, num_intervals ->
        [ meta + [ num_intervals:num_intervals ], cram, crai, intervals, intervals_index ]
    }

// Run process on each interval
PROCESS(cram_intervals, fasta, fasta_fai)

// Gather: Branch by whether intervals were used
vcf_out = PROCESS.out.vcf.branch{
    intervals:    it[0].num_intervals > 1
    no_intervals: it[0].num_intervals <= 1
}

// Merge interval results
vcf_to_merge = vcf_out.intervals
    .map{ meta, vcf -> [ groupKey(meta, meta.num_intervals), vcf ]}
    .groupTuple()

MERGE_VCFS(vcf_to_merge, dict)

// Combine merged and non-interval results, clean up meta
vcf_final = Channel.empty()
    .mix(MERGE_VCFS.out.vcf, vcf_out.no_intervals)
    .map{ meta, vcf -> [ meta - meta.subMap('num_intervals') + [ variantcaller:'toolname' ], vcf ] }
```

---

## Configuration

### Module Configuration Files

Module behavior is controlled via `conf/modules/<tool>.config`:

```groovy
process {
    withName: 'NEWTOOL_PROCESS' {
        ext.args   = { params.newtool_args ?: '' }
        ext.prefix = { "${meta.id}.newtool" }
        publishDir = [
            mode: params.publish_dir_mode,
            path: { "${params.outdir}/variant_calling/newtool/${meta.id}/" },
            pattern: "*{vcf.gz,vcf.gz.tbi}"
        ]
    }
}
```

> **Note:** Older config files wrap process blocks in `if (params.tools && params.tools.split(',').contains('tool'))` guards. Do **not** use this pattern in new code — control which processes run via channel operations (`filter`, `branch`) in the workflow/subworkflow instead. When touching existing config files, remove these guards.

### Resource Labels

Use standard nf-core labels in `conf/base.config`:

| Label                 | CPUs | Memory | Time |
| --------------------- | ---- | ------ | ---- |
| `process_single`      | 1    | 6 GB   | 8 h  |
| `process_low`         | 2    | 12 GB  | 8 h  |
| `process_medium`      | 6    | 36 GB  | 16 h |
| `process_high`        | 12   | 72 GB  | 32 h |
| `process_long`        | -    | -      | 40 h |
| `process_high_memory` | -    | 200 GB | -    |

### Adding New Parameters

1. **Add to `nextflow.config`** with default value:

   ```groovy
   params {
       new_param = false
   }
   ```

2. **Update schema** using nf-core tools:

   ```bash
   nf-core pipelines schema build
   ```

3. **Add validation** if needed in the workflow

---

## Testing

### Test Framework

Sarek uses **nf-test** for testing. Tests are in `tests/` directory.

### Running Tests

```bash
# Run all tests
nf-test test --profile debug,test,docker --verbose

# Run specific test
nf-test test tests/variant_calling_haplotypecaller.nf.test --profile debug,test,docker

# Run with stub mode (faster, no actual execution)
nf-test test tests/default.nf.test --profile debug,test,docker -stub

# Update snapshots when outputs legitimately change
nf-test test tests/my_test.nf.test --profile debug,test,docker --update-snapshot
```

### Test Structure with UTILS.groovy

Tests use a scenario-based pattern:

```groovy
def projectDir = new File('.').absolutePath

nextflow_pipeline {
    name "Test pipeline"
    script "../main.nf"
    tag "pipeline"
    tag "pipeline_sarek"

    def test_scenario = [
        [
            name: "Test scenario name",
            params: [
                input: "${projectDir}/tests/csv/3.0/fastq_single.csv",
                tools: 'haplotypecaller'
            ]
        ],
        [
            name: "Test with stub",
            params: [],
            stub: true
        ],
        [
            name: "Fails with invalid input",
            params: [
                input: "${projectDir}/tests/csv/3.0/vcf_single.csv",
                step: 'annotate',
                vep_cache_version: 1,
                build_only_index: true,
                tools: 'vep'
            ],
            failure: true,
            stdout: "Expected error message"
        ]
    ]

    test_scenario.each { scenario ->
        test(scenario.name, UTILS.get_test(scenario))
    }
}
```

### Test Scenario Options

| Option                         | Description                               |
| ------------------------------ | ----------------------------------------- |
| `name`                         | Test name (descriptive)                   |
| `params`                       | Map of parameters to set                  |
| `stub`                         | Run in stub mode (boolean)                |
| `failure`                      | Expect test to fail (boolean)             |
| `stdout`                       | Expected stdout content for failure tests |
| `gpu`                          | GPU test (adds gpu tag)                   |
| `no_conda`                     | Incompatible with conda                   |
| `include_muse_txt`             | Include MuSE txt in assertions            |
| `include_freebayes_unfiltered` | Include freebayes unfiltered VCFs         |
| `no_vcf_md5sum`                | Use VCF summary instead of md5            |

---

## Documentation

### Documentation Requirements

**Any change that affects pipeline output or adds new functionality MUST include documentation updates.**

### Documentation Files

| File             | Purpose                  | When to Update                                  |
| ---------------- | ------------------------ | ----------------------------------------------- |
| `README.md`      | Pipeline overview        | **New tools** (add to overview/tool list)       |
| `docs/usage.md`  | Usage instructions       | New parameters, new tools, input format changes |
| `docs/output.md` | Output file descriptions | **Any change to outputs**, new tools            |
| `CHANGELOG.md`   | Version history          | Every PR                                        |
| `CITATIONS.md`   | Tool citations           | New tools                                       |
| `docs/images/`   | Metro maps, diagrams     | **New tools**, workflow changes                 |

### New Tool Documentation Checklist

When adding a new tool, you **MUST** update ALL of the following:

1. **`README.md`** - Add tool to the pipeline overview/feature list
2. **`docs/usage.md`** - Document all new parameters and usage instructions
3. **`docs/output.md`** - Document all output files produced by the tool
4. **`docs/images/sarek_subway.*`** - Add tool to the metro map (SVG and PNG)
5. **`CITATIONS.md`** - Add tool citation
6. **`CHANGELOG.md`** - Document the addition

### Output Changes Documentation

Any PR that changes pipeline outputs (new files, changed file names, different content) **MUST** update:

1. **`docs/output.md`** - Reflect the new/changed outputs
2. **`CHANGELOG.md`** - Note the change under appropriate section

### CHANGELOG Format

Follow [Keep a Changelog](https://keepachangelog.com/) format.

**Important conventions:**

- Entries reference the **PR number**, not the issue number:
  ```
  - [#XXX](https://github.com/nf-core/sarek/pull/XXX) - Description of change
  ```
- Use `XXX` as placeholder when no PR exists yet
- The issue number goes in the **PR description body** (for auto-close), not the changelog
- Entries within each section are in **ascending order** by PR number

```markdown
## [Unreleased]

### Added

- [#PR_NUMBER](https://github.com/nf-core/sarek/pull/PR_NUMBER) - Description

### Changed

### Fixed

### Removed

### Dependencies

| Dependency | Old version | New version |
| ---------- | ----------- | ----------- |
| tool_name  | 1.0.0       | 1.1.0       |

### Parameters

| Params        | status |
| ------------- | ------ |
| `--new_param` | New    |

### Developer section

#### Added

#### Changed

#### Fixed

#### Removed
```

### Output Documentation

In `docs/output.md`, document each tool's outputs:

```markdown
### Tool Name

<details markdown="1">
<summary>Output files</summary>

- `path/to/output/`
  - `*.extension`: Description of the file

</details>

Brief description of what this tool produces.
```

---

## Metro Map Updates

### Metro Map Files

Located in `docs/images/`:

- `sarek_subway.svg` / `sarek_subway.png` - Main pipeline flow
- `sarek_indices_subway.svg` / `sarek_indices_subway.png` - Index building flow

### When to Update

- Adding new tools or variant callers
- Adding new preprocessing steps
- Changing the pipeline flow
- Adding new post-processing options

### Update Process

1. Edit the SVG file (use Inkscape or similar)
2. Export to PNG
3. Follow nf-core [design guidelines](https://nf-co.re/developers/design_guidelines)
4. After release, checkout figures from `master` to `dev`:
   ```bash
   git checkout upstream/master -- docs/images/sarek_subway.svg
   git checkout upstream/master -- docs/images/sarek_subway.png
   ```

---

## PR Checklist

### Before Submitting

- [ ] PR targets `dev` branch (not `master`)
- [ ] Code follows Harshil alignment style
- [ ] **Any touched code updated to strict syntax** (explicit closure params, underscore for unused vars)
- [ ] **No new `ext.when` usage** - use channel operations instead
- [ ] **Prefer channel operations** (`filter`, `branch`) over `if` statements for flow control
- [ ] All tests pass: `nf-test test --profile debug,test,docker`
- [ ] Linting passes: `nf-core pipelines lint`
- [ ] No debug mode warnings

### For New Tools

**Code:**

- [ ] Module added/imported correctly
- [ ] Configuration in `conf/modules/<tool>.config`
- [ ] Test added in `tests/`
- [ ] MultiQC config updated (`assets/multiqc_config.yml`) if tool has MultiQC module

**Documentation (ALL required):**

- [ ] `README.md` - Tool added to pipeline overview/feature list
- [ ] `docs/usage.md` - All parameters documented with usage instructions
- [ ] `docs/output.md` - All output files documented
- [ ] `docs/images/sarek_subway.svg` - Tool added to metro map
- [ ] `docs/images/sarek_subway.png` - Exported PNG of updated metro map
- [ ] `CITATIONS.md` - Tool citation added
- [ ] `CHANGELOG.md` - Addition documented

### For New Variant Callers

Adding a variant caller touches **6 locations** — missing any of them causes silent bugs. All of the above "New Tools" items apply, plus:

- [ ] **`nextflow_schema.json`** - Add to the `tools` parameter regex pattern
- [ ] **Dispatcher subworkflow** - Add call and wire outputs into `vcf_all.mix(...)`:
  - `subworkflows/local/bam_variant_calling_germline_all/main.nf` (germline callers)
  - `subworkflows/local/bam_variant_calling_somatic_all/main.nf` (somatic pair callers)
  - `subworkflows/local/bam_variant_calling_tumor_only_all/main.nf` (tumor-only callers)
  - Use channel operations (`filter`, `branch`) to control execution — not `if` blocks (existing `if` blocks are legacy)
- [ ] **`subworkflows/local/post_variantcalling/main.nf`** - Add to `small_variantcallers` list (for SNV callers eligible for normalization/filtering/consensus) or `excluded_variantcallers` (for SV callers). **Forgetting this silently excludes the caller from post-processing.**
- [ ] **Individual subworkflow** - Set `variantcaller` in meta map (e.g., `meta + [variantcaller: 'toolname']`)

### For New Parameters

- [ ] Default value in `nextflow.config`
- [ ] Schema updated: `nf-core pipelines schema build`
- [ ] Validation added (if needed)
- [ ] Documentation in `docs/usage.md`
- [ ] CHANGELOG updated

### For Changes Affecting Pipeline Output

Any PR that changes output files (new files, renamed files, changed content):

- [ ] `docs/output.md` updated to reflect changes
- [ ] `CHANGELOG.md` updated
- [ ] If significant workflow change: metro map updated (`docs/images/sarek_subway.*`)

---

## Common Gotchas

### 1. Forgetting `failOnDuplicate`/`failOnMismatch` on Joins

**Problem:** Silent data loss or incorrect pairing
**Solution:** Always use `join(..., failOnDuplicate: true, failOnMismatch: true)`

### 2. Strelka Produces Two VCFs

**Problem:** Strelka outputs SNV and Indel VCFs separately
**Solution:** Use `groupTuple(size: 2)` then `BCFTOOLS_CONCAT` before downstream processing

### 3. Blocking GroupTuple

**Problem:** `groupTuple()` without size blocks pipeline
**Solution:** Use `groupKey(meta, meta.num_intervals)` when size is known

### 4. Meta Fields Persisting

**Problem:** Temporary meta fields (like `num_intervals`) persist in output
**Solution:** Clean up with `meta - meta.subMap('field_name')` before emit

### 5. DeepVariant Conda

**Problem:** DeepVariant doesn't support Conda
**Solution:** Note in module: `// FIXME Conda is not supported at the moment`

### 6. BWA Memory Requirements

**Problem:** Unexpected OOM errors
**Solution:** BWA requires ~5.37N memory, BWA-MEM2 requires ~28N GB where N = reference size

### 7. Using Deprecated `ext.when` Pattern

**Problem:** Old code uses `ext.when` in config to control module execution
**Solution:** When touching this code, refactor to use channel operations (`filter`, `branch`) to control dataflow. Nextflow is a dataflow language - let the data drive execution. Avoid both `ext.when` AND `if` statements where possible.

### 8. Forgetting to Register a New Variant Caller

**Problem:** New variant caller runs and produces VCFs, but is silently excluded from normalization, filtering, and consensus calling
**Solution:** Must update all 6 registration points — see [For New Variant Callers](#for-new-variant-callers) checklist. The most commonly missed is `post_variantcalling/main.nf`'s `small_variantcallers` list.

### 9. Implicit Variables in Closures

**Problem:** Using implicit `it` makes code harder to read and review
**Solution:** Always use explicit named parameters in closures: `.map { meta, vcf -> ... }` not `.map { it[0], it[1] -> ... }`

---

## Quick Reference

### Essential Commands

```bash
# Run tests
nf-test test --profile debug,test,docker

# Lint pipeline
nf-core pipelines lint

# Update schema
nf-core pipelines schema build

# Install/update module
nf-core modules install <tool>/<subcommand>
nf-core modules update <tool>/<subcommand>
```

### Key Files for Common Changes

| Change Type   | Primary Files                                   |
| ------------- | ----------------------------------------------- |
| New parameter | `nextflow.config`, `nextflow_schema.json`       |
| New tool      | `conf/modules/<tool>.config`, subworkflow, test |
| Bug fix       | Relevant module/subworkflow, test               |
| Documentation | `docs/usage.md`, `docs/output.md`               |
| Any change    | `CHANGELOG.md`                                  |

---

## Getting Help

- **Slack:** [#sarek channel](https://nfcore.slack.com/channels/sarek)
- **GitHub Issues:** [nf-core/sarek/issues](https://github.com/nf-core/sarek/issues)
- **Documentation:** [nf-co.re/sarek](https://nf-co.re/sarek)
