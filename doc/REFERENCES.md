# Download all the Reference files

We are using the [GATK bundle](https://software.broadinstitute.org/gatk/download/bundle) **b37** for most of our reference files.

Following files need to be downloaded (and unzipped):

- '1000G_phase1.indels.b37.vcf'
- '1000G_phase1.indels.b37.vcf.idx'
- 'dbsnp_138.b37.vcf'
- 'dbsnp_138.b37.vcf.idx'
- 'human_g1k_v37_decoy.dict.gz'
- 'human_g1k_v37_decoy.fasta.fai.gz'
- 'human_g1k_v37_decoy.fasta.gz'
- 'Mills_and_1000G_gold_standard.indels.b37.vcf'
- 'Mills_and_1000G_gold_standard.indels.b37.vcf.idx'

This file is on our repo so it will be easy ;-)

- '[centromeres.list](https://raw.githubusercontent.com/SciLifeLab/CAW/master/repeats/centromeres.list)'

The last files are made when indexing human_g1k_v37_decoy.fasta with bwa.

- 'human_g1k_v37_decoy.fasta.pac'
- 'human_g1k_v37_decoy.fasta.amb'
- 'human_g1k_v37_decoy.fasta.ann'
- 'human_g1k_v37_decoy.fasta.bwt'
- 'human_g1k_v37_decoy.fasta.sa'

You can do that by using your own bwa:

```
bwa index human_g1k_v37_decoy.fasta
```

Or using a Docker image containing it:

```
docker run -v `pwd`:/tmp -w /tmp maxulysse/mapreads:1.0 bwa index human_g1k_v37_decoy.fasta
```

The rest of the references files are stored in in [export.uppmax.uu.se](https://export.uppmax.uu.se/b2015110/caw-references/b37/) and also on the repository [CAW-References](https://github.com/MaxUlysse/CAW-References) using [GIT-LFS](https://git-lfs.github.com/) and :

- '1000G_phase3_20130502_SNP_maf0.3.loci'
- 'b37_cosmic_v74.noCHR.sort.4.1.vcf.idx'
- 'b37_cosmic_v74.noCHR.sort.4.1.vcf'

--------------------------------------------------------------------------------

[![](images/SciLifeLab_logo.png "SciLifeLab")][scilifelab-link] [![](images/NGI-final-small.png "NGI")][ngi-link]

[ngi-link]: https://ngisweden.scilifelab.se/
[scilifelab-link]: http://www.scilifelab.se/
