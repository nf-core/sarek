# Use cases

The workflow has two pre-processing options: `mapping` and `recalibrate`.
Using the `mapping` directive one will have a pair of mapped, deduplicated and recalibrated BAM files in the `Preprocessing/Recalibrated/` directory.
This is the usual option you have to give when you are starting from raw FASTQ data:

```bash
nextflow run nf-core/sarek --input mysample.tsv --tools <tool>
```

`mapping` will start by default, you do not have to give any additional parameters, only the TSV file describing the sample (see below).

During the execution of the workflow a `execution_trace.txt`, a `execution_timeline.html` and a `execution_report.html` files are generated automatically.
These files contain statistics about resources used and processes finished.
If you start a new workflow or restart/resume a sample, the previous version will be renamed as `execution_trace.txt.1`, `execution_timeline.html.1` and `execution_report.html.1` respectively.
Also, older version are renamed with incremented numbers.

## Starting from raw FASTQ - pair of FASTQ files

The workflow should be started in this case with the smallest set of options as written above:

```bash
nextflow run nf-core/sarek --input mysample.tsv --tools <tool>
```

The TSV file should look like:

```text
SUBJECT_ID  XX    0    SAMPLE_ID    1    /samples/normal_1.fastq.gz    /samples/normal_2.fastq.gz
```

See the [input files documentation](input.md) for more information.

## Starting from raw FASTQ - a directory with normal sample only

The `--input` option can be also used to point Sarek to a directory with FASTQ files:

```bash
nextflow run nf-core/sarek --input path/to/FASTQ/files --tools <tool>
```

The given directory is searched recursively for FASTQ files that are named `*_R1_*.fastq.gz`, and a matching pair with the same name except `_R2_` instead of `_R1_` is expected to exist alongside.
All of the found FASTQ files are considered to belong to the sample.
Each FASTQ file pair gets its own read group (`@RG`) in the resulting BAM file.

### Metadata when using `--input` with a directory

When using `--input` with a directory, the metadata about the sample that are written to the BAM header in the `@RG` tag are determined in the following way.

- The sample name (`SM`) is derived from the the last component of the path given to `--input`.
That is, you should make sure that that directory has a meaningful name! For example, with `--input=/my/fastqs/sample123`, the sample name will be `sample123`.
- The read group id is set to *flowcell.samplename.lane*.
The flowcell id and lane number are auto-detected from the name of the first read in the FASTQ file.

## Starting from raw FASTQ - pair of FASTQ files for tumor/normal samples

The workflow command line is just the same as before, but the TSV contains extra lines.
You can see the second column is used to distinguish normal and tumor samples.
You can add as many relapse samples as many you have, providing their name in the third column is different.
Each will be compared to the normal one-by-one.
Usually there are more read groups - sequencing lanes - for a single sequencing run, and in a flowcell different lanes have to be recalibrated separately.
This is captured in the TSV file only in the following manner, adding read group numbers or IDs in the fourth column.
All lanes belonging to the same Sample will be merged together after the FASTQ pairs are mapped to the reference genome.
Obviously, if you do not have relapse samples, you can leave out the two last lines.

```text
SUBJECT_ID  XX    0    SAMPLE_ID_N    1    /samples/normal1_1.fastq.gz    /samples/normal1_2.fastq.gz
SUBJECT_ID  XX    0    SAMPLE_ID_N    2    /samples/normal2_1.fastq.gz    /samples/normal2_2.fastq.gz
SUBJECT_ID  XX    1    SAMPLE_ID_T    3    /samples/tumor3_1.fastq.gz    /samples/tumor3_2.fastq.gz
SUBJECT_ID  XX    1    SAMPLE_ID_T    4    /samples/tumor4_1.fastq.gz    /samples/tumor4_2.fastq.gz
SUBJECT_ID  XX    1    SAMPLE_ID_T    5    /samples/tumor5_1.fastq.gz    /samples/tumor5_2.fastq.gz
SUBJECT_ID  XX    1    SAMPLE_ID_R    7    /samples/relapse7_1.fastq.gz    /samples/relapse7_2.fastq.gz
SUBJECT_ID  XX    1    SAMPLE_ID_R    9    /samples/relapse9_1.fastq.gz    /samples/relapse9_2.fastq.gz
```

See the [input files documentation](input.md) for more information.

## Starting from recalibration

```bash
nextflow run nf-core/sarek --input mysample.tsv --step recalibrate --tools <tool>
```

And the corresponding TSV file should be like:
Obviously, if you do not have tumor or relapse samples, you can leave out the two last lines.

```text
SUBJECT_ID  XX    0    SAMPLE_ID_N    /samples/SAMPLE_ID_N.bam    /samples/SAMPLE_ID_N.bai /samples/SAMPLE_ID_N.recal.table
SUBJECT_ID  XX    1    SAMPLE_ID_T    /samples/SAMPLE_ID_T.bam    /samples/SAMPLE_ID_T.bai /samples/SAMPLE_ID_T.recal.table
SUBJECT_ID  XX    1    SAMPLE_ID_R    /samples/SAMPLE_ID_R.bam    /samples/SAMPLE_ID_R.bai /samples/SAMPLE_ID_R.recal.table
```

See the [input files documentation](input.md) for more information.

## Starting from a recalibrated BAM file

At this step we are assuming that all the required preprocessing is over, we only want to run variant callers or other tools using recalibrated BAM files.

```bash
nextflow run nf-core/sarek --step variantcalling --tools <tool>
```

And the corresponding TSV file should be like:

```text
SUBJECT_ID  XX    0    SAMPLE_ID_N    /samples/SAMPLE_ID_N.bam    /samples/SAMPLE_ID_N.bai
SUBJECT_ID  XX    1    SAMPLE_ID_T    /samples/SAMPLE_ID_T.bam    /samples/SAMPLE_ID_T.bai
SUBJECT_ID  XX    1    SAMPLE_ID_R    /samples/SAMPLE_ID_R.bam    /samples/SAMPLE_ID_R.bai
```

See the [input files documentation](input.md) for more information.

If you want to restart a previous run of the pipeline, you may not have a recalibrated BAM file.
In this case, you need to start with `--step=recalibrate` (see previous section).

## Using Sarek with targeted (whole exome or panel) sequencing data

The recommended flow for targeted sequencing data is to use the workflow as it is, but also provide a BED file containing targets for all steps using the `--targetBED` option.
The workflow will pick up these intervals, and activate the `--exome` flag to process deeper coverage.
It is adviced to pad the variant calling regions (exons or the target) to some extent before submitting to the workflow.
To add the target BED file configure the flow like:

```bash
nextflow run nf-core/sarek --tools haplotypecaller,strelka,mutect2 --targetBED targets.bed --input my_panel.tsv
```
