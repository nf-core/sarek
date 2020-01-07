# Genomes and reference files

## AWS iGenomes

Sarek is using [AWS iGenomes](https://ewels.github.io/AWS-iGenomes/), which facilitate storing and sharing references.
Sarek currently uses `GRCh38` by default.
`GRCh37`,  `GRCh38` and `GRCm38` are available with `--genome GRCh37`, `--genome GRCh38` or `--genome GRCm38` respectively with any profile using the `conf/igenomes.config` file, or you can specify it with `-c conf/igenomes.config`.
Use `--genome smallGRCh37` to map against a small reference genome based on GRCh37.
Settings in `igenomes.config` can be tailored to your needs.

### Intervals

To speed up some preprocessing and variant calling processes, the reference is chopped into smaller pieces.
The intervals are chromosomes cut at their centromeres (so each chromosome arm processed separately) also additional unassigned contigs.
We are ignoring the `hs37d5` contig that contains concatenated decoy sequences.
Parts of preprocessing and variant calling are done by these intervals, and the different resulting files are then merged.
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

If no intervals files are specified, one will be automatically generated following:

```bash
awk -v FS='\t' -v OFS='\t' '{ print \$1, \"0\", \$2 }' <REFERENCE>.fasta.fai > <REFERENCE>.bed
```

To disable this feature, please use [`--no_intervals`](usage.md#--no_intervals)

### Working with whole exome (WES) or panel data

The `--targetBED` parameter does _not_  imply that the workflow is running alignment or variant calling only for the supplied targets.
Instead, we are aligning for the whole genome, and selecting variants only at the very end by intersecting with the provided target file.
Adding every exon as an interval in case of WES can generate >200K processes or jobs, much more forks, and similar number of directories in the Nextflow work directory.
Furthermore, primers and/or baits are not 100% specific, (certainly not for MHC and KIR, etc.), quite likely there going to be reads mapping to multiple locations.
If you are certain that the target is unique for your genome (all the reads will certainly map to only one location), and aligning to the whole genome is an overkill, better to change the reference itself.

### Working with other genomes

> :warning: This is a new feature, in active development, so usage could change.

Sarek can also do limited preprocessing from any genome, providing a `fasta` file as a reference genome, followed by limited variant calling using `mpileup`, `Manta` and `Strelka`.

Limited support for `TAIR10`, `EB2`, `UMD3.1`, `bosTau8`, `WBcel235`, `ce10`, `CanFam3.1`, `canFam3`, `GRCz10`, `danRer10`, `BDGP6`, `dm6`, `EquCab2`, `equCab2`, `EB1`, `Galgal4`, `galGal4`, `Gm01`, `hg38`, `hg19`, `Mmul_1`, `mm10`, `IRGSP-1.0`, `CHIMP2.1.4`, `panTro4`, `Rnor_6.0`, `rn6`, `R64-1-1`, `sacCer3`, `EF2`, `Sbi1`, `Sscrofa10.2`, `susScr3`, `AGPv3`.
