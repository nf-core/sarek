# Genomes and reference files

Sarek currently uses GRCh38 by default.
The settings are in `genomes.config`, they can be tailored to your needs.
The [`buildReferences.nf`](#buildreferencesnf) script is used to build the indexes for the reference test.

## GRCh37

Use `--genome GRCh37` to map against GRCh37.
Before doing so and if you are not on UPPMAX, you need to adjust the settings in `genomes.config` to your needs.

### GATK bundle

To get the needed files, download the [GATK bundle for GRCh37](ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/).

The following files need to be downloaded:

- 242c0df2a698a76fc43bdd938ba57c62 - '1000G\_phase1.indels.b37.vcf.gz'
- 00b0e74e4a13536dd6c0728c66db43f3 - 'dbsnp\_138.b37.vcf.gz'
- dd05833f18c22cc501e3e31406d140b0 - 'human\_g1k\_v37\_decoy.fasta.gz'
- a0764a80311aee369375c5c7dda7e266 - 'Mills\_and\_1000G\_gold\_standard.indels.b37.vcf.gz'

### Other files for GRCh37

From our repo, get the [`intervals` list file](https://raw.githubusercontent.com/SciLifeLab/Sarek/master/repeats/wgs_calling_regions.grch37.list).
More information about this file in the [intervals documentation](INTERVALS.md)

Description of how to generate the Loci file used in the ASCAT process is described [here](https://github.com/SciLifeLab/Sarek/blob/master/docs/ASCAT.md).

## GRCh38

Use `--genome GRCh38` to map against GRCh38.
Before doing so and if you are not on UPPMAX, you need to adjust the settings in `genomes.config` to your needs.

To get the needed files, download the GATK bundle for GRCh38 from [ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/](ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/).
You can also download the required files from the Google Cloud mirror link [here](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0).

The MD5SUM of `Homo_sapiens_assembly38.fasta` included in that file is 7ff134953dcca8c8997453bbb80b6b5e.

If you download the data from the FTP servers `beta/` directory, which seems to be an older version of the bundle, only `Homo_sapiens_assembly38.known_indels.vcf` is needed.
Also, you can omit `dbsnp_138_` and `dbsnp_144` files as we use `dbsnp_146`.
The old ones also use the wrong chromosome naming convention.
The Google Cloud mirror has all data in the `v0` directory, but requires you to remove the `resources_broad_hg38_v0_` prefixes from all files.

The following files need to be downloaded:

- 3884c62eb0e53fa92459ed9bff133ae6 - 'Homo_sapiens_assembly38.dict'
- 7ff134953dcca8c8997453bbb80b6b5e - 'Homo_sapiens_assembly38.fasta'
- b07e65aa4425bc365141756f5c98328c - 'Homo_sapiens_assembly38.fasta.64.alt'
- e4dc4fdb7358198e0847106599520aa9 - 'Homo_sapiens_assembly38.fasta.64.amb'
- af611ed0bb9487fb1ba4aa1a7e7ad21c - 'Homo_sapiens_assembly38.fasta.64.ann'
- d41d8cd98f00b204e9800998ecf8427e - 'Homo_sapiens_assembly38.fasta.64.bwt'
- 178862a79b043a2f974ef10e3877ef86 - 'Homo_sapiens_assembly38.fasta.64.pac'
- 91a5d5ed3986db8a74782e5f4519eb5f - 'Homo_sapiens_assembly38.fasta.64.sa'
- f76371b113734a56cde236bc0372de0a - 'Homo_sapiens_assembly38.fasta.fai'
- 14cc588a271951ac1806f9be895fb51f - 'Homo_sapiens_assembly38.known_indels.vcf.gz'
- 1a55fdfa6533ae5cbc70e8188e779229 - 'Homo_sapiens_assembly38.known_indels.vcf.gz.tbi'
- 2e02696032dcfe95ff0324f4a13508e3 - 'Mills_and_1000G_gold_standard.indels.hg38.vcf.gz'
- 4c807e2cbe0752c0c44ac82ff3b52025 - 'Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi'

If you just downloaded the `Homo_sapiens_assembly38.fasta.gz` file, you would need to do:

```
gunzip Homo_sapiens_assembly38.fasta.gz
bwa index -6 Homo_sapiens_assembly38.fasta
```

Description of how to generate the Loci file used in the ASCAT process is described [here](https://github.com/SciLifeLab/Sarek/blob/master/docs/ASCAT.md).

## smallGRCh37

Use `--genome smallGRCh37` to map against a small reference genome based on GRCh37.
`smallGRCh37` is the default genome for the testing profile (`-profile testing`).

## AWS iGenomes
Sarek is using [AWS iGenomes](https://ewels.github.io/AWS-iGenomes/), which facilitate storing and sharing references.
Both `GRCh37` and `GRCh38` are available with `--genome GRCh37` or `--genome GRCh38` respectively with any profile using the `conf/igenomes.config` file (eg.: `awsbatch`, or `btb`), or you can specify it with `-c conf/igenomes.config`, it contains all data previously detailed.

## buildReferences.nf

The `buildReferences.nf` script can download and build the files needed for smallGRCh37, or build the references for GRCh37/smallGRCh37.

### `--refDir`

Use `--refDir <path to references>` to specify where are the files to process.

```
nextflow run buildReferences.nf --refDir <path to references> --genome <genome>
```

### `--genome`

Same parameter used for other scripts.

- GRCh37
- GRCh38 (not yet available)
- smallGRCh37
