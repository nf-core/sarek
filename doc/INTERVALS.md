# Intervals

To speed up the Variant Calling processes, the reference is chopped into smaller pieces. The Variant Calling is done by this intervals, and the different resulting VCFs are then merged. This can parallelize the Variant Calling processes, and push down the variant calling wall clock time significanlty.

## GRCh37

The file `wgs_calling_regions.grch37.list` contains the intervals which are chopped up at the centromeres.

## GRCh38

The file `wgs_calling_regions_sorted.hg38.list` contains the intervals which are a list of non-N regions, sorted by their size, to enhance wall clock time.

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
