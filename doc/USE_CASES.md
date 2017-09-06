# Use cases

The workflow has four pre-processing options: `mapping`, `realign`, `recalibrate` and `variantcalling`. Using the `mapping` directive one will have a pair of mapped, deduplicated and recalibrated BAM files in the `Preprocessing/Recalibrated/` directory. Furthermore, during this process a deduplicated BAM file is created in the `Preprocessing/NonRealigned/` directory. This is the usual option you have to give when you are starting from raw FASTQ data:

```bash
nextflow run SciLifeLab/CAW --sample mysample.tsv
```

`mapping` will start by default, you do not have to give any additional parameters, only the TSV file describing the sample (see below).

In the [genomes.config](../configuration/genomes.config) configuration file we are defining the intervals file as well, this is used to define regions for variant call and realignment (in a scatter and gather fashion when possible). The intervals are chromosomes cut at their centromeres (so each chromosome arm processed separately) also additional unassigned contigs. We are ignoring the hs37d5 contig that contains concatenated decoy sequences.

During the execution of the workflow a `trace.txt` and a `timeline.html` file is generated automatically. These files contain statistics about resources used and processes finished. If you start a new flow or restart/resume a sample, the previous version will be renamed as `trace.txt.1` and `timeline.html.1` respectively. Also, older version are renamed with incremented numbers.

## Starting from raw FASTQ - having pair of FASTQ files for normal and tumor samples (one lane for each sample)

The workflow should be started in this case with the smallest set of options as written above:

```bash
nextflow run SciLifeLab/CAW --sample mysample.tsv
```

The TSV file should have at least two tab-separated lines:

```
SUBJECT_ID  XX    0    SAMPLEIDN    1    /samples/normal_1.fastq.gz    /samples/normal_2.fastq.gz
SUBJECT_ID  XX    1    SAMPLEIDT    1    /samples/tumor_1.fastq.gz    /samples/tumor_2.fastq.gz
```

The columns are:

1. Subject id
2. status: 0 if normal, 1 if tumor
3. sample id: actual text representation of the type of the sample
4. read group ID: it is irrelevant in this simple case, should to be 1
5. first set of reads
6. second set of reads

## Starting from raw FASTQ - having pair of FASTQ files for normal, tumor and relapse samples (one lane for each sample)

The workflow command line is just the same as before, but the TSV contains an extra line. You can see the second column is used to distinguish normal and tumor, there is no extra identifier for relapse. You can add as many relapse samples as many you have, providing their name in the third column is different. Each will be compared to the normal one-by-one.

```
SUBJECT_ID  XX    0    SAMPLEIDN    1    /samples/normal_1.fastq.gz    /samples/normal_2.fastq.gz
SUBJECT_ID  XX    1    SAMPLEIDT    1    /samples/tumor_1.fastq.gz    /samples/tumor_2.fastq.gz
SUBJECT_ID  XX    1    SAMPLEIDR    1    /samples/relapse_1.fastq.gz    /samples/relapse_2.fastq.gz
```

## Starting from raw FASTQ - having multiple lanes (reads groups)

Usually there are more read groups - sequencing lanes - for a single sequencing run, and in a flowcell different lanes have to be recalibrated separately. This is captured in the TSV file only in the following manner, adding read group numbers or IDs in the fourth column. Obviously, if you do not have relapse samples, you can leave out those lines.

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
nextflow run SciLifeLab/CAW --sample mysample.tsv --step realign
```

And the corresponding TSV file should be like:

```
SUBJECT_ID  XX    0    SAMPLEIDN    /samples/SAMPLEIDN.bam    /samples/SAMPLEIDN.bai
SUBJECT_ID  XX    1    SAMPLEIDT    /samples/SAMPLEIDT.bam    /samples/SAMPLEIDT.bai
```

At the end of this step you should have recalibrated BAM files in the `Preprocessing/Recalibrated/` directory.

## Starting from recalibration

If the BAM files were realigned together, you can start from recalibration:

```bash
nextflow run SciLifeLab/CAW --sample mysample.tsv --step recalibrate
```

And the corresponding TSV file should be like:

```
SUBJECT_ID  XX    0    SAMPLEIDN    /samples/SAMPLEIDN.bam    /samples/SAMPLEIDN.bai /samples/SAMPLEIDN.recal.table
SUBJECT_ID  XX    1    SAMPLEIDT    /samples/SAMPLEIDT.bam    /samples/SAMPLEIDT.bai /samples/SAMPLEIDT.recal.table
```

## Starting from a recalibrated BAM file

At this step we are assuming that all the required preprocessing (alignment, deduplication, ..., recalibration) is over, we only want to run variant callers or other tools using recalibrated BAMs.

```bash
nextflow run SciLifeLab/CAW --sample mysample.tsv --step variantcalling --tools <tool>
```

And the corresponding TSV file should be like:

```
SUBJECT_ID  XX    0    SAMPLEIDN    /samples/SAMPLEIDN.bam    /samples/SAMPLEIDN.bai
SUBJECT_ID  XX    1    SAMPLEIDT    /samples/SAMPLEIDT.bam    /samples/SAMPLEIDT.bai
```

If you want to restart a previous run of the pipeline, you may not have a recalibrated BAM file. This is the case if HaplotypeCaller was the only tool (recalibration is done on-the-fly with HaplotypeCaller to improve performance and save space). In this case, you need to start with `--step=recalibrate` (see previous section).

--------------------------------------------------------------------------------

[![](images/SciLifeLab_logo.png "SciLifeLab")][scilifelab-link]
[![](images/NGI_logo.png "NGI")][ngi-link]
[![](images/NBIS_logo.png "NBIS")][nbis-link]

[nbis-link]: https://www.nbis.se/
[ngi-link]: https://ngisweden.scilifelab.se/
[scilifelab-link]: https://www.scilifelab.se/
