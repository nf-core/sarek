# Genomes and reference files

Sarek currently uses GRCh38 by default. The settings are in `genomes.config`, they can be tailored to your needs. The [`buildReferences.nf`](#buildreferencesnf) script can be use to build the indexes based on the reference files.

## GRCh37

Use `--genome GRCh37` to map against GRCh37. Before doing so and if you are not on UPPMAX, you need to adjust the settings in `genomes.config` to your needs.

### GATK bundle

To get the needed files, download the [GATK bundle for GRCh37](ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/).

The following files need to be downloaded:

- 242c0df2a698a76fc43bdd938ba57c62 - '1000G\_phase1.indels.b37.vcf.gz'
- 00b0e74e4a13536dd6c0728c66db43f3 - 'dbsnp\_138.b37.vcf.gz'
- dd05833f18c22cc501e3e31406d140b0 - 'human\_g1k\_v37\_decoy.fasta.gz'
- a0764a80311aee369375c5c7dda7e266 - 'Mills\_and\_1000G\_gold\_standard.indels.b37.vcf.gz'

### Other files

From our repo, get the [`intervals` list file](https://raw.githubusercontent.com/SciLifeLab/Sarek/master/repeats/wgs_calling_regions.grch37.list). More information about this file in the [intervals documentation](INTERVALS.md)

You can create your own cosmic reference for any human reference as specified below.

### COSMIC files

To annotate with COSMIC variants during MuTect1/2 Variant Calling you need to create a compatible VCF file.
Download the coding and non-coding VCF files from [COSMIC](http://cancer.sanger.ac.uk/cosmic/download) and
process them with the [Create\_Cosmic.sh](https://github.com/SciLifeLab/Sarek/tree/master/scripts/Create_Cosmic.sh)
script. The script requires a fasta index `.fai`, of the reference file you are using.

Example:

```bash
samtools faidx human_g1k_v37_decoy.fasta
sh Create_Cosmic.sh human_g1k_v37_decoy.fasta.fai
```

Note: CosmicCodingMuts.vcf.gz & CosmicNonCodingVariants.vcf.gz must be in same folder as Create\_Cosmic.sh when executed.

To index the resulting VCF file use [igvtools](https://software.broadinstitute.org/software/igv/igvtools).

```bash
igvtools index <cosmicvxx.vcf>
```

## GRCh38

Use `--genome GRCh38` to map against GRCh38. Before doing so and if you are not on UPPMAX, you need to adjust the settings in `genomes.config` to your needs.

To get the needed files, download the GATK bundle for GRCh38 from [ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/](ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/).

The MD5SUM of `Homo_sapiens_assembly38.fasta` included in that file is 7ff134953dcca8c8997453bbb80b6b5e.

From the `beta/` directory, which seems to be an older version of the bundle, only `Homo_sapiens_assembly38.known_indels.vcf` is needed. Also, you can omit `dbsnp_138_` and `dbsnp_144` files as we use `dbsnp_146`. The old ones also use the wrong chromosome naming convention.

Afterwards, the following needs to be done:

```
gunzip Homo_sapiens_assembly38.fasta.gz
bwa index -6 Homo_sapiens_assembly38.fasta
```

## smallGRCh37

Use `--genome smallGRCh37` to map against a small reference genome based on GRCh37. `smallGRCh37` is the default genome for the testing profile (`-profile testing`).

## buildReferences.nf

The `buildReferences.nf` script can download and build the files needed for smallGRCh37, or build the references for GRCh37/smallGRCh37.

### `--download`

Only with `--genome smallGRCh37`. If this option is specify, the [`smallRef`](https://github.com/szilvajuhos/smallRef) repository will be automatically downloaded from GitHub. Not to be used on UPPMAX cluster Bianca or on similarly secured clusters where such things are not working/allowed.

```
nextflow run buildReferences.nf --download --genome smallGRCh37
```

### `--refDir`

Use `--refDir <path to smallRef>` to specify where are the files to process.

```
nextflow run buildReferences.nf --refDir <path to smallRef> --genome <genome>
```

### `--genome`

Same parameter used for other scripts.

- GRCh37
- GRCh38 (not yet available)
- smallGRCh37

--------------------------------------------------------------------------------

[![](images/SciLifeLab_logo.png "SciLifeLab")][scilifelab-link]
[![](images/NGI_logo.png "NGI")][ngi-link]
[![](images/NBIS_logo.png "NBIS")][nbis-link]

[nbis-link]: https://www.nbis.se/
[ngi-link]: https://ngisweden.scilifelab.se/
[scilifelab-link]: https://www.scilifelab.se/
