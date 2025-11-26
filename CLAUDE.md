# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

nf-core/sarek is a Nextflow DSL2 pipeline for detecting germline and somatic variants from whole genome or targeted sequencing data. It supports multiple aligners (bwa-mem, bwa-mem2, dragmap, parabricks) and numerous variant callers (GATK HaplotypeCaller/Mutect2, Strelka, Manta, DeepVariant, Freebayes, Sentieon, etc.).

## Build and Test Commands

### Running the pipeline locally

```bash
# Basic test run
nextflow run main.nf -profile test,docker --outdir results

# Run with debug profile for process selector warnings
nextflow run main.nf -profile debug,test,docker --outdir results
```

### Running nf-test

```bash
# Run all tests
nf-test test --profile debug,test,docker --verbose

# Run a specific test file
nf-test test tests/default.nf.test --profile debug,test,docker --verbose

# Run tests matching a tag
nf-test test --tag cpu --profile test,docker

# Run stub tests (faster, no real execution)
nf-test test --profile test,docker -stub
```

### Linting

```bash
# Run nf-core linting
nf-core pipelines lint

# Update parameter schema after adding new params
nf-core pipelines schema build
```

## Architecture

### Entry Points

- `main.nf` - Main entry point that:
  - Loads genome parameters via `getGenomeAttribute()`
  - Calls `NFCORE_SAREK` workflow which orchestrates preparation and main analysis
  - Handles `PREPARE_GENOME`, `PREPARE_INTERVALS`, and annotation cache setup
  - Passes prepared references to the `SAREK` workflow

### Core Workflow

- `workflows/sarek/main.nf` - Main analysis workflow containing:
  - Input validation and routing based on `step` parameter (mapping, markduplicates, prepare_recalibration, recalibrate, variant_calling, annotate)
  - Preprocessing: `FASTQ_PREPROCESS_GATK` or `FASTQ_PREPROCESS_PARABRICKS`
  - Variant calling: `BAM_VARIANT_CALLING_GERMLINE_ALL`, `BAM_VARIANT_CALLING_TUMOR_ONLY_ALL`, `BAM_VARIANT_CALLING_SOMATIC_ALL`
  - Post-processing: `POST_VARIANTCALLING` (filtering, normalization, consensus calling)
  - Annotation: `VCF_ANNOTATE_ALL` (snpEff, VEP)

### Subworkflows Organization (`subworkflows/local/`)

- **Preprocessing**: `bam_markduplicates/`, `bam_applybqsr/`, `bam_baserecalibrator/`
- **Variant Calling**: Named by pattern `bam_variant_calling_{germline|somatic|tumor_only}_{tool}/`
- **CSV Generation**: `channel_*_create_csv/` - Create restart samplesheets at checkpoints
- **Post-processing**: `vcf_normalization/`, `vcf_consensus/`, `post_variantcalling/`

### Configuration Structure (`conf/`)

- `base.config` - Default resource allocations with process labels
- `test.config` - Minimal test dataset configuration
- `modules/*.config` - Tool-specific process configurations (publishDir, ext.args, etc.)
- `igenomes.config` - Reference genome paths for common genomes

### Test Structure

- Tests use nf-test framework with plugins: `nft-bam`, `nft-utils`, `nft-vcf`
- `tests/lib/UTILS.groovy` - Shared test helper providing `get_test()` and `get_assertion()` functions
- Test scenarios defined as maps with params, tags (cpu/gpu, conda compatibility), stub options
- Assertions check: versions YAML, stable filenames, BAM/CRAM reads MD5, VCF variants MD5

## Key Parameters

- `--step` - Pipeline entry point: mapping (default), markduplicates, prepare_recalibration, recalibrate, variant_calling, annotate
- `--tools` - Comma-separated variant callers/annotators to run
- `--skip_tools` - Tools to skip (e.g., markduplicates, baserecalibrator)
- `--aligner` - bwa-mem (default), bwa-mem2, dragmap, parabricks
- `--wes` - Set true for exome/targeted sequencing data

## Input Format

Samplesheet CSV with columns: `patient,sample,lane,fastq_1,fastq_2` (for FASTQs) or with BAM/CRAM inputs depending on start step.

## Module Updates

Modules from nf-core/modules live in `modules/nf-core/`. Update with:

```bash
nf-core modules update <module_name>
```

Local modules are in `modules/local/`.
