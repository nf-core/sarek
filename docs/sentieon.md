# nf-core/sarek: Usage with sentieon <!-- omit in toc -->

- [Introduction](#introduction)
- [Sentieon Analysis Pipelines & Tools](#sentieon-analysis-pipelines--tools)
  - [Alignment](#alignment)
  - [Germline SNV/INDEL Variant Calling - DNAseq](#germline-snvindel-variant-calling---dnaseq)
  - [Germline SNV/INDEL Variant Calling - DNAscope](#germline-snvindel-variant-calling---dnascope)
  - [Somatic SNV/INDEL Variant Calling - TNscope](#somatic-snvindel-variant-calling---tnscope)
  - [Structural Variant Calling](#structural-variant-calling)
- [usage](#usage)
  - [--sentieon](#--sentieon)
  - [--tools](#--tools)

## Introduction

[Sentieon](https://www.sentieon.com/) is a commercial solution to process genomics data with high computing efficiency, fast turnaround time, exceptional accuracy, and 100% consistency.

If [Sentieon](https://www.sentieon.com/) is available, use this `--sentieon` params to enable with Sarek to use some Sentieon Analysis Pipelines & Tools.

Please refer to the [nf-core/configs](https://github.com/nf-core/configs#adding-a-new-pipeline-specific-config) repository on how to make a pipeline-specific configuration file based on the [munin-sarek specific configuration file](https://github.com/nf-core/configs/blob/master/conf/pipeline/sarek/munin.config).

Or ask us on the [nf-core Slack](http://nf-co.re/join/slack) on the following channels: [#sarek](https://nfcore.slack.com/channels/sarek) or [#configs](https://nfcore.slack.com/channels/configs).

## Sentieon Analysis Pipelines & Tools

The following Sentieon Analysis Pipelines & Tools are available within Sarek.

### Alignment

> Sentieon BWA matches BWA-MEM with > 2X speedup.

This tool is enabled by default within Sarek if `--sentieon` is specified and if the pipeline is started with the `mapping` [step](usage.md#--step).

### Germline SNV/INDEL Variant Calling - DNAseq

> Precision FDA award-winning software.
> Matches GATK 3.3-4.1, and without down-sampling.
> Results up to 10x faster and 100% consistent every time.

This tool is enabled within Sarek if `--sentieon` is specified and if `--tools DNAseq` is specified cf [--tools](#--tools).

### Germline SNV/INDEL Variant Calling - DNAscope

> Improved accuracy and genome characterization.
> Machine learning enhanced filtering producing top variant calling accuracy.

This tool is enabled within Sarek if `--sentieon` is specified and if `--tools DNAscope` is specified cf [--tools](#--tools).

### Somatic SNV/INDEL Variant Calling - TNscope

> Winner of ICGC-TCGA DREAM challenge.
> Improved accuracy, machine learning enhanced filtering.
> Supports molecular barcodes and unique molecular identifiers.

This tool is enabled within Sarek if `--sentieon` is specified and if `--tools TNscope` is specified cf [--tools](#--tools).

### Structural Variant Calling

> Germline and somatic SV calling, including translocations, inversions, duplications and large INDELs

This tool is enabled within Sarek if `--sentieon` is specified and if `--tools DNAscope` is specified cf [--tools](#--tools).

## usage

### --sentieon

Adds the following tools for the [`--tools`](#--tools) options: `DNAseq`, `DNAscope` and `TNscope`.

### --tools

For main usage of tools, follow the [usage/tools](usage.md#--tools) documentation.

With `--sentieon` the following tools options are also available within Sarek: `DNAseq`, `DNAscope` and `TNscope`.
