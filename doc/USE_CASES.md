# Use cases

The workflow has three pre-processing options: `mapping`, `realign` and `recalibrate`. Using the `mapping` directive one will have a pair of mapped, deduplicated and recalibrated BAM files in the `Preprocessing/Recalibrated/` directory. Furthermore, during this process a deduplicated BAM file is created in the `Preprocessing/NonRealigned/` directory. This is the usual option you have to give when you are starting from raw FASTQ data:

```bash
nextflow run SciLifeLab/Sarek/main.nf --sample mysample.tsv
nextflow run SciLifeLab/Sarek/germlineVC.nf --tools <tool>
nextflow run SciLifeLab/Sarek/somaticVC.nf --tools <tool> # For somatic only
nextflow run SciLifeLab/Sarek/annotate.nf --tool <tool> --annotateVCF myfile.vcf # For somatic only
nextflow run SciLifeLab/Sarek/runMultiQC.nf
```

`mapping` will start by default, you do not have to give any additional parameters, only the TSV file describing the sample (see below).

In the [genomes.config](https://raw.githubusercontent.com/SciLifeLab/Sarek/master/configuration/genomes.config) configuration file we are defining the intervals file as well, this is used to define regions for variant call and realignment (in a scatter and gather fashion when possible). The intervals are chromosomes cut at their centromeres (so each chromosome arm processed separately) also additional unassigned contigs. We are ignoring the hs37d5 contig that contains concatenated decoy sequences.

During the execution of the workflow a `Sarek-trace.txt`, a `Sarek-timeline.html` and a `Sarek-report.html` files are generated automatically. These files contain statistics about resources used and processes finished. If you start a new workflow or restart/resume a sample, the previous version will be renamed as `Sarek-trace.txt.1`, `Sarek-timeline.html.1` and `Sarek-report.html.1` respectively. Also, older version are renamed with incremented numbers.

## Starting from raw FASTQ - pair of FASTQ files

The workflow should be started in this case with the smallest set of options as written above:

```bash
nextflow run SciLifeLab/Sarek/main.nf --sample mysample.tsv
nextflow run SciLifeLab/Sarek/germlineVC.nf --tools <tool>
nextflow run SciLifeLab/Sarek/runMultiQC.nf
```

The TSV file should have at least one tab-separated lines:

```
SUBJECT_ID  XX    0    SAMPLEID    1    /samples/normal_1.fastq.gz    /samples/normal_2.fastq.gz
```

The columns are:

1. Subject id
2. status: 0 if normal, 1 if tumor
3. sample id: actual text representation of the type of the sample
4. read group ID: it is irrelevant in this simple case, should to be 1
5. first set of reads
6. second set of reads

## Starting from raw FASTQ on a normal sample only (with `--sampleDir`)

The `--sampleDir` option can be used to point Sarek to a directory with FASTQ files:
```bash
nextflow run SciLifeLab/Sarek/main.nf --sampleDir path/to/FASTQ/files
nextflow run SciLifeLab/Sarek/germlineVC.nf --tools <tool>
nextflow run SciLifeLab/Sarek/runMultiQC.nf
```
The given directory is searched recursively for FASTQ files that are named `*_R1_*.fastq.gz`, and a matching pair with the same name except `_R2_` instead of `_R1_` is expected to exist alongside. All of the found FASTQ files are considered to belong to the sample. Each FASTQ file pair gets its own read group (`@RG`) in the resulting BAM file.

### Metadata when using `--sampleDir`

When using `--sampleDir`, the metadata about the sample that are written to the BAM header in the `@RG` tag are determined in the following way.

- The sample name (`SM`) is derived from the the last component of the path given to `--sampleDir`. That is, you should make sure that that directory has a meaningful name! For example, with `--sampleDir=/my/fastqs/sample123`, the sample name will be `sample123`.
- The read group id is set to *flowcell.samplename.lane*. The flowcell id and lane number are auto-detected from the name of the first read in the FASTQ file.

## Starting from raw FASTQ - having pair of FASTQ files for tumor/normal samples (one lane for each sample)

The workflow command line is just the same as before, but the TSV contains extra lines. You can see the second column is used to distinguish normal and tumor samples. You can add as many relapse samples as many you have, providing their name in the third column is different. Each will be compared to the normal one-by-one. Obviously, if you do not have relapse samples, you can leave out this last line.

```
SUBJECT_ID  XX    0    SAMPLEIDN    1    /samples/normal_1.fastq.gz    /samples/normal_2.fastq.gz
SUBJECT_ID  XX    1    SAMPLEIDT    1    /samples/tumor_1.fastq.gz    /samples/tumor_2.fastq.gz
SUBJECT_ID  XX    1    SAMPLEIDR    1    /samples/relapse_1.fastq.gz    /samples/relapse_2.fastq.gz
```

## Starting from raw FASTQ - having multiple lanes (reads groups)

Usually there are more read groups - sequencing lanes - for a single sequencing run, and in a flowcell different lanes have to be recalibrated separately. This is captured in the TSV file only in the following manner, adding read group numbers or IDs in the fourth column.

```
SUBJECT_ID  XX    0    SAMPLEID    1    /samples/normal1_1.fastq.gz    /samples/normal1_2.fastq.gz
SUBJECT_ID  XX    0    SAMPLEID    2    /samples/normal2_1.fastq.gz    /samples/normal2_2.fastq.gz
```

## Starting from raw FASTQ - having multiple lanes (reads groups) for tumor/normal samples

Usually there are more read groups - sequencing lanes - for a single sequencing run, and in a flowcell different lanes have to be recalibrated separately. This is captured in the TSV file only in the following manner, adding read group numbers or IDs in the fourth column. Obviously, if you do not have relapse samples, you can leave out those last two lines.

```
SUBJECT_ID  XX    0    SAMPLEIDN    1    /samples/normal1_1.fastq.gz    /samples/normal1_2.fastq.gz
SUBJECT_ID  XX    0    SAMPLEIDN    2    /samples/normal2_1.fastq.gz    /samples/normal2_2.fastq.gz
SUBJECT_ID  XX    1    SAMPLEIDT    3    /samples/tumor3_1.fastq.gz    /samples/tumor3_2.fastq.gz
SUBJECT_ID  XX    1    SAMPLEIDT    4    /samples/tumor4_1.fastq.gz    /samples/tumor4_2.fastq.gz
SUBJECT_ID  XX    1    SAMPLEIDT    5    /samples/tumor5_1.fastq.gz    /samples/tumor5_2.fastq.gz
SUBJECT_ID  XX    1    SAMPLEIDR    7    /samples/relapse7_1.fastq.gz    /samples/relapse7_2.fastq.gz
SUBJECT_ID  XX    1    SAMPLEIDR    9    /samples/relapse9_1.fastq.gz    /samples/relapse9_2.fastq.gz
```

## Starting from realignement

NGI Production in the previous years delivered many preprocessed samples; these BAM files are not recalibrated. To have BAMs suitable for variant calling, realignement of pairs is necessary:

```bash
nextflow run SciLifeLab/Sarek/main.nf --sample mysample.tsv --step realign
nextflow run SciLifeLab/Sarek/germlineVC.nf --tools <tool>
nextflow run SciLifeLab/Sarek/runMultiQC.nf
```

And the corresponding TSV file should be like:

```
SUBJECT_ID  XX    0    SAMPLEID    /samples/SAMPLEIDN.bam    /samples/SAMPLEIDN.bai
```

At the end of this step you should have recalibrated BAM files in the `Preprocessing/Recalibrated/` directory.

## Starting from realignement for tumor/normal samples

NGI Production in the previous years delivered many preprocessed samples; these BAM files are not recalibrated. To have BAMs suitable for variant calling, realignement of pairs is necessary:

```bash
nextflow run SciLifeLab/Sarek/main.nf --sample mysample.tsv --step realign
nextflow run SciLifeLab/Sarek/germlineVC.nf --tools <tool>
nextflow run SciLifeLab/Sarek/somaticVC.nf --tools <tool>
nextflow run SciLifeLab/Sarek/annotate.nf --tool <tool> --annotateVCF myfile.vcf
nextflow run SciLifeLab/Sarek/runMultiQC.nf

```

And the corresponding TSV file should be like (obviously, if you do not have relapse samples, you can leave out this last line):

```
SUBJECT_ID  XX    0    SAMPLEIDN    /samples/SAMPLEIDN.bam    /samples/SAMPLEIDN.bai
SUBJECT_ID  XX    1    SAMPLEIDT    /samples/SAMPLEIDT.bam    /samples/SAMPLEIDT.bai
SUBJECT_ID  XX    1    SAMPLEIDR    /samples/SAMPLEIDT.bam    /samples/SAMPLEIDR.bai
```

At the end of this step you should have recalibrated BAM files in the `Preprocessing/Recalibrated/` directory.

## Starting from recalibration for tumor/normal samples

If the BAM files were realigned together, you can start from recalibration:

```bash
nextflow run SciLifeLab/Sarek/main.nf --sample mysample.tsv --step recalibrate
nextflow run SciLifeLab/Sarek/germlineVC.nf --tools <tool>
nextflow run SciLifeLab/Sarek/somaticVC.nf --tools <tool>
nextflow run SciLifeLab/Sarek/annotate.nf --tool <tool> --annotateVCF myfile.vcf
nextflow run SciLifeLab/Sarek/runMultiQC.nf
```

And the corresponding TSV file should be like (obviously, if you do not have relapse samples, you can leave out this last line):

```
SUBJECT_ID  XX    0    SAMPLEIDN    /samples/SAMPLEIDN.bam    /samples/SAMPLEIDN.bai /samples/SAMPLEIDN.recal.table
SUBJECT_ID  XX    1    SAMPLEIDT    /samples/SAMPLEIDT.bam    /samples/SAMPLEIDT.bai /samples/SAMPLEIDT.recal.table
SUBJECT_ID  XX    1    SAMPLEIDR    /samples/SAMPLEIDR.bam    /samples/SAMPLEIDR.bai /samples/SAMPLEIDR.recal.table
```

## Starting from a recalibrated BAM file

At this step we are assuming that all the required preprocessing is over, we only want to run variant callers or other tools using recalibrated BAMs.

```bash
nextflow run SciLifeLab/Sarek/germlineVC.nf --tools <tool>
nextflow run SciLifeLab/Sarek/runMultiQC.nf
```

And the corresponding TSV file should be like:

```
SUBJECT_ID  XX    0    SAMPLEIDN    /samples/SAMPLEIDN.bam    /samples/SAMPLEIDN.bai
```

If you want to restart a previous run of the pipeline, you may not have a recalibrated BAM file. This is the case if HaplotypeCaller was the only tool (recalibration is done on-the-fly with HaplotypeCaller to improve performance and save space). In this case, you need to start with `--step=recalibrate` (see previous section).

## Starting from a recalibrated BAM file for tumor/normal samples

At this step we are assuming that all the required preprocessing is over, we only want to run variant callers or other tools using recalibrated BAMs.

```bash
nextflow run SciLifeLab/Sarek/germlineVC.nf --tools <tool>
nextflow run SciLifeLab/Sarek/somaticVC.nf --tools <tool>
nextflow run SciLifeLab/Sarek/annotate.nf --tool <tool> --annotateVCF myfile.vcf
nextflow run SciLifeLab/Sarek/runMultiQC.nf
```

And the corresponding TSV file should be like (obviously, if you do not have relapse samples, you can leave out this last line):

```
SUBJECT_ID  XX    0    SAMPLEIDN    /samples/SAMPLEIDN.bam    /samples/SAMPLEIDN.bai
SUBJECT_ID  XX    1    SAMPLEIDT    /samples/SAMPLEIDT.bam    /samples/SAMPLEIDT.bai
SUBJECT_ID  XX    1    SAMPLEIDR    /samples/SAMPLEIDR.bam    /samples/SAMPLEIDR.bai
```

If you want to restart a previous run of the pipeline, you may not have a recalibrated BAM file. This is the case if HaplotypeCaller was the only tool (recalibration is done on-the-fly with HaplotypeCaller to improve performance and save space). In this case, you need to start with `--step=recalibrate` (see previous section).

--------------------------------------------------------------------------------

[![](images/SciLifeLab_logo.png "SciLifeLab")][scilifelab-link]
[![](images/NGI_logo.png "NGI")][ngi-link]
[![](images/NBIS_logo.png "NBIS")][nbis-link]

[nbis-link]: https://www.nbis.se/
[ngi-link]: https://ngisweden.scilifelab.se/
[scilifelab-link]: https://www.scilifelab.se/
