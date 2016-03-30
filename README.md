# CAW
Cancer Analysis Workflow Prototype at SciLifeLab

## Quick start

Install [nextflow][nextflow] and then run

```
nextflow somaticPipe.nf
```

which will run the full workflow on a small test data set (se below).

To run you own sample, override with:

```
nextflow somaticPipe.nf --sample <SAMPLE>
```

Note that the workflow assumes that there are input files for <SAMPLE>:
./data/<SAMPLE>.tumor_R1.fastq.gz
./data/<SAMPLE>.tumor_R2.fastq.gz
./data/<SAMPLE>.normal_R1.fastq.gz
./data/<SAMPLE>.normal_R2.fastq.gz


## Test data

For quick, easy testing of the workflow, we need a small, but representative
data set as input. We went for a full depth coverage, small chunk of the genome.
"Seba" suggested using [TCGA benchmarking data][TCGA], which is already downloaded
to Uppmax.

Grabbing 100kb from chr1, both tumor and normal, sorting by read name:
```
$ samtools view -bh /sw/data/uppnex/ToolBox/TCGAbenchmark/f0eaa94b-f622-49b9-8eac-e4eac6762598/G15511.HCC1143_BL.1.bam 1:100000-200000 | samtools sort -n -m 4G - tcga.cl.normal.part.nsort.bam
$ samtools view -bh /sw/data/uppnex/ToolBox/TCGAbenchmark/ad3d4757-f358-40a3-9d92-742463a95e88/G15511.HCC1143.1.bam 1:100000-200000 | samtools sort -n -m 4G - tcga.cl.tumor.part.nsort.bam
```

Then convert to fastq:
```
$ bedtools bamtofastq -i tcga.cl.normal.part.nsort.bam -fq tcga.cl.normal.part.nsort_R1.fastq -fq2 tcga.cl.normal.part.nsort_R2.fastq
$ bedtools bamtofastq -i tcga.cl.tumor.part.nsort.bam -fq tcga.cl.tumor.part.nsort_R1.fastq -fq2 tcga.cl.tumor.part.nsort_R2.fastq
```
And finally gzipping...


[TCGA]: https://cghub.ucsc.edu/datasets/benchmark_download.html
[nextflow]: http://www.nextflow.io
