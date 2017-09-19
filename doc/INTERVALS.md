# Intervals

To speed up the Variant Calling processes, the reference is chopped into smaller pieces. The Variant Calling is done by this intervals, and the different resulting VCFs are then merged. This can parallelize the Variant Calling processes, and push down the variant calling wall clock time significanlty.

The calling intervals can be defined using a `.list` or a `.bed` file. A `.list` file contains one interval per line in the format `chromosome:start-end` (1-based coordinates).

When the intervals file is in BED format, the file must be a tab-separated text file with one interval per line. There must be at least three columns: chromosome, start, and end. In BED format, the coordinates are 0-based, so the interval `chrom:1-10` becomes `chrom<tab>0<tab>10`.

Additionally, the "score" column of the BED file can be used to provide an estimate of how many seconds it will take to call variants on that interval. The fourth column remains unused. Example (the fields would actually be tab-separated, this is not shown here):

    chr1  10000  207666 NA  47.3

This indicates that variant calling on the interval chr1:10001-207666 takes approximately 47.3 seconds.

The runtime estimate is used in two different ways. First, when there are multiple consecutive intervals in the file that take little time to compute, they are processed as a single job, thus reducing the number of processes that needs to be spawned. Second, the jobs with largest processing time are started first, which reduces wall-clock time. If no runtime is given, a time of 1000 nucleotides per second is assumed. Actual figures vary from 2 nucleotides/second to 30000 nucleotides/second.

## GRCh37

The file `wgs_calling_regions.grch37.list` contains the intervals which are chopped up at the centromeres.

## GRCh38

The file `wgs_calling_regions_sorted.hg38.list` contains the intervals which are a list of non-N regions, sorted by their size, to enhance wall clock time.

The file `wgs_calling_regions.hg38.bed` contains the intervals which are a list of non-N regions in a BED format, with the fifth column containing runtime estimate used to enhance wall clock time.

## smallGRCh37

The file `tiny_GRCh37.list` is based on the `GRCh37` intervals list file but tailored for the `smallReferences`.

## smallGRCh38

The file `tiny_GRCh38.list` contains the intervals from the `tiny_GRCh37.list` lifted over from `GRCh37` to `GRCh38` in the following way:

```bash
module load bioinfo-tools samtools/1.4 bwa
samtools faidx /sw/data/uppnex/ToolBox/ReferenceAssemblies/hg38make/bundle/2.8/b37/human_g1k_v37_decoy.fasta $(< tiny.list) > intervals.fasta
bwa mem /sw/data/uppnex/ToolBox/hg38bundle/Homo_sapiens_assembly38.fasta intervals.fasta > intervals.sam
grep -v '^@' intervals.sam | awk '{printf("%s:%d-%d\n", $3, $4, $4+$6-1)}' > tiny-GRCh38.list
```

--------------------------------------------------------------------------------

[![](images/SciLifeLab_logo.png "SciLifeLab")][scilifelab-link]
[![](images/NGI_logo.png "NGI")][ngi-link]
[![](images/NBIS_logo.png "NBIS")][nbis-link]

[nbis-link]: https://www.nbis.se/
[ngi-link]: https://ngisweden.scilifelab.se/
[scilifelab-link]: https://www.scilifelab.se/
