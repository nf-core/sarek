# Genomes and reference files

## AWS iGenomes

Sarek is using [AWS iGenomes](https://ewels.github.io/AWS-iGenomes/), which facilitate storing and sharing references.
Sarek currently uses `GRCh38` by default.
Both `GRCh37` and `GRCh38` are available with `--genome GRCh37` or `--genome GRCh38` respectively with any profile using the `conf/igenomes.config` file, or you can specify it with `-c conf/igenomes.config`.
Use `--genome smallGRCh37` to map against a small reference genome based on GRCh37.
Settings in `igenomes.config` can be tailored to your needs.

### Intervals

To speed up some preprocessing and variant calling processes, the reference is chopped into smaller pieces.
The intervals are chromosomes cut at their centromeres (so each chromosome arm processed separately) also additional unassigned contigs.
We are ignoring the hs37d5 contig that contains concatenated decoy sequences.
Parts of preprocessing and variant calling are done by this intervals, and the different resulting files are then merged.
This can parallelize processes, and push down wall clock time significantly.

The calling intervals can be defined using a `.list` or a `.bed` file.
A `.list` file contains one interval per line in the format `chromosome:start-end` (1-based coordinates).

When the intervals file is in BED format, the file must be a tab-separated text file with one interval per line.
There must be at least three columns: chromosome, start, and end.
In BED format, the coordinates are 0-based, so the interval `chrom:1-10` becomes `chrom<tab>0<tab>10`.

Additionally, the "score" column of the BED file can be used to provide an estimate of how many seconds it will take to call variants on that interval.
The fourth column remains unused.
Example (the fields would actually be tab-separated, this is not shown here):

`chr1  10000  207666 NA  47.3`

This indicates that variant calling on the interval chr1:10001-207666 takes approximately 47.3 seconds.

The runtime estimate is used in two different ways.
First, when there are multiple consecutive intervals in the file that take little time to compute, they are processed as a single job, thus reducing the number of processes that needs to be spawned.
Second, the jobs with largest processing time are started first, which reduces wall-clock time.
If no runtime is given, a time of 1000 nucleotides per second is assumed.
Actual figures vary from 2 nucleotides/second to 30000 nucleotides/second.

## build.nf

The [`build.nf`](#buildnf) script is used to build reference needed for smallGRCh37.

```bash
nextflow run build.nf
```
