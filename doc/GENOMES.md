# Genomes

CAW currently uses GRCh37 by default. Support for GRCh38 is not fully working!
The settings are in `genomes.config`, they can be tailored to your needs.

## GRCh37

...


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

### buildReferences.nf
The `buildReferences.nf` script can build/dowload the files needed for testing.
The files are located within the [`smallRef` repository](https://github.com/szilvajuhos/smallRef).

If the `--dowload` option is specify, the `smallRef` repo will be automatically downloaded from github. Do not use on UPPMAX cluster Bianca where such things are not allowed.
```
nextflow run buildReferences.nf --download --genome smallGRCh37
```
You can use instead `--refDir <path to smallRef>` to process
```
nextflow run buildReferences.nf --refDir <path to smallRef> --genome smallGRCh37
```
No need to specify the `--genome` option if `smallGRCh37` is already specified in your config files or in your profile.
