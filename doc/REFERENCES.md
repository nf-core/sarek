# Genomes and specific reference files

CAW currently uses GRCh37 by default. Support for GRCh38 is not fully working yet!
The settings are in `genomes.config`, they can be tailored to your needs. The [`buildReferences.nf`](#buildReferences.nf) script can be use to build the indexes based on the reference files.

## GRCg37
Use `--genome GRCh37` to map against GRCh37. Before doing so and if you
are not on Uppmax, you need to adjust the settings in `genomes.config` to your
needs.

### GATK bundle
To get the needed files, download the [GATK bundle for GRCh37](ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/).

The following files need to be downloaded:

- '1000G_phase1.indels.b37.vcf.gz' - 242c0df2a698a76fc43bdd938ba57c62
- 'dbsnp_138.b37.vcf.gz' - 00b0e74e4a13536dd6c0728c66db43f3
- 'human_g1k_v37_decoy.fasta.gz' - dd05833f18c22cc501e3e31406d140b0
- 'Mills_and_1000G_gold_standard.indels.b37.vcf.gz' - a0764a80311aee369375c5c7dda7e266

### COSMIC files

To annotate with COSMIC variants during mutect1/2 variant calling you need to create a compatible VCF file.
Download the coding and non-coding VCF files from [COSMIC](http://cancer.sanger.ac.uk/cosmic/download) and process them with the [Create_Cosmic.sh](https://github.com/SciLifeLab/CAW/tree/master/scripts/Create_Cosmic.sh) script. The script requires a fasta index, .fai, of the reference file you are using.

Example:

```
samtools faidx human_g1k_v37_decoy.fasta
sh Create_Cosmic.sh human_g1k_v37_decoy.fasta.fai
```

Note: CosmicCodingMuts.vcf.gz & CosmicNonCodingVariants.vcf.gz must be in same folder as Create_Cosmic.sh when executed.

To index the resulting VCF file use [igvtools](https://software.broadinstitute.org/software/igv/igvtools).

```
igvtools index COSMICv##.vcf
```

### Other files
From our repo, get the '[centromeres.list](https://raw.githubusercontent.com/SciLifeLab/CAW/master/repeats/centromeres.list)' file.

The rest of the references files are stored in in [export.uppmax.uu.se](https://export.uppmax.uu.se/b2015110/caw-references/b37/) and also on the repository [CAW-References](https://github.com/MaxUlysse/CAW-References) using [GIT-LFS](https://git-lfs.github.com/):

- '1000G_phase3_20130502_SNP_maf0.3.loci'
- 'b37_cosmic_v74.noCHR.sort.4.1.vcf'

You can create your own cosmic reference for any human reference as specified below.

## GRCh38

Use `--genome=GRCh38` to map against GRCh38. Before doing so and if you
are not on Uppmax, you need to adjust the settings in `genomes.config` to your
needs.

To get the needed files, download the GATK bundle for GRCh38 from
<ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/>.

The MD5SUM of `Homo_sapiens_assembly38.fasta` included in that file is
7ff134953dcca8c8997453bbb80b6b5e.

From the `beta/` directory, which seems to be an older version of the bundle,
only Homo_sapiens_assembly38.known_indels.vcf* is needed. Also, you can omit
dbsnp_138* and dbsnp_144 files as we use dbsnp_146. The old ones also use the
wrong chromosome naming convention.

Afterwards, the following needs to be done:

    gunzip Homo_sapiens_assembly38.fasta.gz
    bwa index -6 Homo_sapiens_assembly38.fasta
    awk '!/^@/{printf("%s:%d-%d\n", $1, $2, $3)}' wgs_calling_regions.hg38.interval_list > wgs_calling_regions.hg38.list

## smallGRCh37
Use `--genome smallGRCh37` to map against a small reference genome based on GRCh37. `smallGRCh37` is the default genome for the testing profile (`-profile testing`).

## buildReferences.nf
The `buildReferences.nf` script can dowload and build the files needed for smallGRCh37, or build the references for GRCh37/smallGRCh37.

### `--download`
Only with `--genome smallGRCh37`
If the `--dowload` option is specify, the [`smallRef` repository](https://github.com/szilvajuhos/smallRef). repo will be automatically downloaded from github. Do not use on UPPMAX cluster Bianca or on similar clusters where such things are not allowed.
```
nextflow run buildReferences.nf --download --genome smallGRCh37
```
### `--refDir`
Use `--refDir <path to smallRef>` to process
```
nextflow run buildReferences.nf --refDir <path to smallRef> --genome <genome>
### `--genome`
Same parameter used for `main.nf`
- GRCh37
- GRch38 (not yet supported)
- smallGRCh37

--------------------------------------------------------------------------------

[![](images/SciLifeLab_logo.png "SciLifeLab")][scilifelab-link] [![](images/NGI-final-small.png "NGI")][ngi-link]

[ngi-link]: https://ngisweden.scilifelab.se/
[scilifelab-link]: http://www.scilifelab.se/
