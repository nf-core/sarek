# Input Documentation

## General information about the TSV files

Input files for Sarek can be specified using a TSV (Tab Separated Values) file given to the `--input` command (Note, the delimiter is the tab (`\t`) character).
There are different kinds of TSV files that can be used as input, depending on the input files available (FASTQ, uBAM, BAM...).
For all possible TSV files, described in the next sections, here is an explanation of what the columns refer to:

- `subject` designates the subject, it should be the ID of the subject, and it must be unique for each subject, but one subject can have multiple samples (e.g.
normal and tumor)
- `sex` are the sex chromosomes of the subject, (XX or XY)
- `status` is the status of the measured sample, (0 for Normal or 1 for Tumor)
- `sample` designates the sample, it should be the ID of the sample (it is possible to have more than one tumor sample for each subject, i.e.
a tumor and a relapse), it must be unique, but samples can have multiple lanes (which will later be merged)
- `lane` is used when the sample is multiplexed on several lanes, it must be unique for each lane in the same sample (but does not need to be the original lane name), and must contain at least one character
- `fastq1` is the path to the first pair of the FASTQ file
- `fastq2` is the path to the second pair of the FASTQ file
- `bam` is the path to the BAM file
- `bai` is the path to the BAM index file
- `recaltable` is the path to the recalibration table
- `mpileup` is the path to the mpileup file

It is recommended to add the absolute path of the files, but relative path should also work.

All examples are given for a normal/tumor pair.
If no tumors are listed in the TSV file, then the workflow will proceed as if it is a normal sample instead of a normal/tumor pair, producing the germline Variant Calling results only.

Sarek will output results in a different directory for each sample.
If multiple samples are specified in the TSV file, Sarek will consider all files to be from different samples.
Multiple TSV files can be specified if the path is enclosed in quotes.

Output from Variant Calling and/or Annotation will be in a specific directory for each sample (or normal/tumor pair if applicable).

## Starting from the mapping step

When starting from the mapping step (`--step mapping`), the first step of Sarek, the input can have three different forms:

- A TSV file containing the sample metadata and the path to the paired-end FASTQ files.
- The path to a directory containing the FASTQ files
- A TSV file containing the sample metadata and the path to the unmapped BAM (uBAM) files.

### Providing a TSV file with the path to FASTQ files

The TSV file to start with the step mapping with paired-end FASTQs should contain the columns:

`subject sex status sample lane fastq1 fastq2`

In this sample for the normal case there are 3 read groups, and 2 for the tumor.

| | | | | | |
|-|-|-|-|-|-|
|G15511|XX|0|C09DFN|C09DF_1|/path/to/C09DFACXX111207.1_1.fastq.gz|/path/to/C09DFACXX111207.1_2.fastq.gz|
|G15511|XX|0|C09DFN|C09DF_2|/path/to/C09DFACXX111207.2_1.fastq.gz|/path/to/C09DFACXX111207.2_2.fastq.gz|
|G15511|XX|0|C09DFN|C09DF_3|/path/to/C09DFACXX111207.3_1.fastq.gz|/path/to/C09DFACXX111207.3_2.fastq.gz|
|G15511|XX|1|D0ENMT|D0ENM_1|/path/to/D0ENMACXX111207.1_1.fastq.gz|/path/to/D0ENMACXX111207.1_2.fastq.gz|
|G15511|XX|1|D0ENMT|D0ENM_2|/path/to/D0ENMACXX111207.2_1.fastq.gz|/path/to/D0ENMACXX111207.2_2.fastq.gz|

### Providing the path to a FASTQ directory

Input files for Sarek can be specified using the path to a FASTQ directory given to the `--input` command only with the `mapping` step.

```bash
nextflow run nf-core/sarek --input </path/To/Directory> ...
```

#### Input FASTQ file name best practices

The input folder, containing the FASTQ files for one subject (ID) should be organized into one sub-folder for every sample.
All FASTQ files for that sample should be collected here.

```text
ID
+--sample1
+------sample1_lib_flowcell-index_lane_R1_1000.fastq.gz
+------sample1_lib_flowcell-index_lane_R2_1000.fastq.gz
+------sample1_lib_flowcell-index_lane_R1_1000.fastq.gz
+------sample1_lib_flowcell-index_lane_R2_1000.fastq.gz
+--sample2
+------sample2_lib_flowcell-index_lane_R1_1000.fastq.gz
+------sample2_lib_flowcell-index_lane_R2_1000.fastq.gz
+--sample3
+------sample3_lib_flowcell-index_lane_R1_1000.fastq.gz
+------sample3_lib_flowcell-index_lane_R2_1000.fastq.gz
+------sample3_lib_flowcell-index_lane_R1_1000.fastq.gz
+------sample3_lib_flowcell-index_lane_R2_1000.fastq.gz
```

FASTQ filename structure:

- `sample_lib_flowcell-index_lane_R1_1000.fastq.gz` and
- `sample_lib_flowcell-index_lane_R2_1000.fastq.gz`

Where:

- `sample` = sample id
- `lib` = identifier of library preparation
- `flowcell` = identifier of flow cell for the sequencing run
- `lane` = identifier of the lane of the sequencing run

Read group information will be parsed from FASTQ file names according to this:

- `RGID` = "sample_lib_flowcell_index_lane"
- `RGPL` = "Illumina"
- `PU` = sample
- `RGLB` = lib

### Providing a TSV file with the paths to uBAM files

The TSV file for starting the mapping from uBAM files should contain the columns:

- `subject sex status sample lane bam`

In this sample for the normal case there are 3 read groups, and 2 for the tumor.

| | | | | | |
|-|-|-|-|-|-|
|G15511|XX|0|C09DFN|C09DF_1|/path/to/C09DFAC_1.bam|
|G15511|XX|0|C09DFN|C09DF_2|/path/to/C09DFAC_2.bam|
|G15511|XX|0|C09DFN|C09DF_3|/path/to/C09DFAC_3.bam|
|G15511|XX|1|D0ENMT|D0ENM_1|/path/to/D0ENMAC_1.bam|
|G15511|XX|1|D0ENMT|D0ENM_2|/path/to/D0ENMAC_2.bam|

## Starting from the BAM prepare recalibration step

To start from the preparation of the recalibration step (`--step prepare_recalibration`), a TSV file for a normal/tumor pair needs to be given as input containing the paths to the non recalibrated but already mapped BAM files.
The TSV needs to contain the following columns:

- `subject sex status sample bam bai`

The same way, if you have non recalibrated BAMs and their indexes, you should use a structure like:

| | | | | | |
|-|-|-|-|-|-|
|G15511|XX|0|C09DFN|/path/to/G15511.C09DFN.md.bam|/path/to/G15511.C09DFN.md.bai|
|G15511|XX|1|D0ENMT|/path/to/G15511.D0ENMT.md.bam|/path/to/G15511.D0ENMT.md.bai|

When starting Sarek from the mapping step, a TSV file is generated automatically after the `MarkDuplicates` process.
This TSV file is stored under `results/Preprocessing/TSV/duplicates_marked_no_table.tsv` and can be used to restart Sarek from the non-recalibrated BAM files.
Using the parameter `--step prepare_recalibration` will automatically take this file as input.

Additionally, individual TSV files for each sample (`duplicates_marked_no_table_[SAMPLE].tsv`) can be found in the same directory.

If `--skip_markduplicates` has been specified, the TSV file for this step will be slightly different:

| | | | | | |
|-|-|-|-|-|-|
|G15511|XX|0|C09DFN|/path/to/G15511.C09DFN.bam|/path/to/G15511.C09DFN.bai|
|G15511|XX|1|D0ENMT|/path/to/G15511.D0ENMT.bam|/path/to/G15511.D0ENMT.bai|

When starting Sarek from the mapping step with `--skip_markduplicates`, a TSV file is generated automatically after the `Mapping` processes.
This TSV file is stored under `results/Preprocessing/TSV/mapped.tsv` and can be used to restart Sarek from the non-recalibrated BAM files.
Using the parameter `--step recalibrate --skip_markduplicates` will automatically take this file as input.

Additionally, individual TSV files for each sample (`mapped_[SAMPLE].tsv`) can be found in the same directory.

## Starting from the BAM recalibration step

To start from the recalibration step (`--step recalibrate`), a TSV file for a normal/tumor pair needs to be given as input containing the paths to the non recalibrated but already mapped BAM files.
The TSV needs to contain the following columns:

- `subject sex status sample bam bai recaltable`

The same way, if you have non recalibrated BAMs, their indexes and their recalibration tables, you should use a structure like:

| | | | | | | |
|-|-|-|-|-|-|-|
|G15511|XX|0|C09DFN|/path/to/G15511.C09DFN.md.bam|/path/to/G15511.C09DFN.md.bai|/path/to/G15511.C09DFN.recal.table|
|G15511|XX|1|D0ENMT|/path/to/G15511.D0ENMT.md.bam|/path/to/G15511.D0ENMT.md.bai|/path/to/G15511.D0ENMT.recal.table|

When starting Sarek from the mapping step, a TSV file is generated automatically after the `BaseRecalibrator` processes.
This TSV file is stored under `results/Preprocessing/TSV/duplicates_marked.tsv` and can be used to restart Sarek from the non-recalibrated BAM files.
Using `--step recalibrate` will automatically take this file as input.

Additionally, individual TSV files for each sample (`duplicates_marked_[SAMPLE].tsv`) can be found in the same directory.

If `--skip_markduplicates` has been specified, the TSV file for this step will be slightly different:

| | | | | | | |
|-|-|-|-|-|-|-|
|G15511|XX|0|C09DFN|/path/to/G15511.C09DFN.bam|/path/to/G15511.C09DFN.bai|/path/to/G15511.C09DFN.recal.table|
|G15511|XX|1|D0ENMT|/path/to/G15511.D0ENMT.bam|/path/to/G15511.D0ENMT.bai|/path/to/G15511.D0ENMT.recal.table|

When starting Sarek from the mapping step with `--skip_markduplicates`, a TSV file is generated automatically after the `BaseRecalibrator` processes.
This TSV file is stored under `results/Preprocessing/TSV/mapped_no_duplicates_marked.tsv` and can be used to restart Sarek from the non-recalibrated BAM files.
Using `--step recalibrate` will automatically take this file as input.

Additionally, individual TSV files for each sample (`mapped_no_duplicates_marked_[SAMPLE].tsv`) can be found in the same directory.

## Starting from the variant calling step

A TSV file for a normal/tumor pair with recalibrated BAM files and their indexes can be provided to start Sarek from the variant calling step (`--step variantcalling`).
The TSV file should contain the columns:

- `subject sex status sample bam bai`

Here is an example for two samples from the same subject:

| | | | | | |
|-|-|-|-|-|-|
|G15511|XX|0|C09DFN|/path/to/G15511.C09DFN.recal.bam|/path/to/G15511.C09DFN.recal.bai|
|G15511|XX|1|D0ENMT|/path/to/G15511.D0ENMT.recal.bam|/path/to/G15511.D0ENMT.recal.bai|

When starting Sarek from the mapping or recalibrate steps, a TSV file is generated automatically after the recalibration processes.
This TSV file is stored under `results/Preprocessing/TSV/recalibrated.tsv` and can be used to restart Sarek from the recalibrated BAM files.
Using the parameter `--step variantcalling` will automatically take this file as input.

Additionally, individual TSV files for each sample (`recalibrated_[SAMPLE].tsv`) can be found in the same directory.

## Starting from the mpileup file with the Control-FREEC step

To start from the Control-FREEC step (`--step Control-FREEC`), a TSV file for a normal/tumor pair needs to be given as input containing the paths to the mpileup files.
The TSV needs to contain the following columns:

- `subject sex status sample mpileup`

Here is an example for one normal/tumor pair from one subjects:

| | | | | |
|-|-|-|-|-|
|G15511|XX|0|C09DFN|/path/to/G15511.C09DFN.pileup|
|G15511|XX|1|D0ENMT|/path/to/G15511.D0ENMT.pileup|

When starting Sarek from the Control-FREEC step, a TSV file is generated automatically after the `mpileup` process.
This TSV file is stored under `results/VariantCalling/TSV/control-freec_mpileup.tsv` and can be used to restart Sarek from the mpileup files.
Using the parameter `--step Control-FREEC` will automatically take this file as input.

Additionally, individual TSV files for each sample (`control-freec_mpileup_[SAMPLE].tsv`) can be found in the same directory.

## VCF files for annotation

Input files for Sarek can be specified using the path to a VCF directory given to the `--input` command only with the annotation step (`--step annotate`).
As Sarek will use `bgzip` and `tabix` to compress and index VCF files annotated, it expects VCF files to be sorted.
Multiple VCF files can be specified, using a [glob path](https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob), if enclosed in quotes.
For example:

```bash
nextflow run nf-core/sarek --step annotate --input "results/VariantCalling/*/{HaplotypeCaller,Manta,Mutect2,Strelka,TIDDIT}/*.vcf.gz" ...
```
