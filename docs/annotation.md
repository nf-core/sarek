# Annotation

## Tools

With Sarek, annotation is done using `snpEff`, `VEP`, or even both consecutively:

- `--tools snpEff`
  - To annotate using `snpEff`
- `--tools VEP`
  - To annotate using `VEP`
- `--tools snpEff,VEP`
  - To annotate using `snpEff` and `VEP`
- `--tools merge`
  - To annotate using `snpEff` followed by `VEP`

VCF produced by Sarek will be annotated if `snpEff` or `VEP` are specified with the `--tools` command.
As Sarek will use `bgzip` and `tabix` to compress and index VCF files annotated, it expects VCF files to be sorted.

In these examples, all command lines will be launched starting with `--step annotate`.
It can of course be started directly from any other step instead.

## Using genome specific containers

Sarek has already designed containers with `snpEff` and `VEP` files for Human (`GRCh37`, `GRCh38`), Mouse (`GRCm38`), Dog (`CanFam3.1`) and Roundworm (`WBcel235`).
Default settings will run using these containers.

The main Sarek container has also `snpEff` and `VEP` installed, but without the cache files that can be downloaded separately.

## Download cache

A Nextflow helper script has been designed to help downloading `snpEff` and `VEP` caches.
Such files are meant to be shared between multiple users, so this script is mainly meant for people administrating servers, clusters and advanced users.

```bash
nextflow run download_cache.nf --snpeff_cache </path/to/snpEff/cache> --snpeff_db <snpEff DB version> --genome <GENOME>
nextflow run download_cache.nf --vep_cache </path/to/VEP/cache> --species <species> --vep_cache_version <VEP cache version> --genome <GENOME>
```

## Using downloaded cache

Both `snpEff` and `VEP` enable usage of cache.
If cache is available on the machine where Sarek is run, it is possible to run annotation using cache.
You need to specify the cache directory using `--snpeff_cache` and `--vep_cache` in the command lines or within configuration files.
The cache will only be used when `--annotation_cache` and cache directories are specified (either in command lines or in a configuration file).

Example:

```bash
nextflow run nf-core/sarek --tools snpEff --step annotate --sample <file.vcf.gz> --snpeff_cache </path/to/snpEff/cache> --annotation_cache
nextflow run nf-core/sarek --tools VEP --step annotate --sample <file.vcf.gz> --vep_cache </path/to/VEP/cache> --annotation_cache
```

## Using VEP CADD plugin

To enable the use of the VEP CADD plugin:

- Download the CADD files
- Specify them (either on the command line, like in the example or in a configuration file)
- use the `--cadd_cache` flag

Example:

```bash
nextflow run nf-core/sarek --step annotate --tools VEP --sample <file.vcf.gz> --cadd_cache \
    --cadd_indels </path/to/CADD/cache/InDels.tsv.gz> \
    --cadd_indels_tbi </path/to/CADD/cache/InDels.tsv.gz.tbi> \
    --cadd_wg_snvs </path/to/CADD/cache/whole_genome_SNVs.tsv.gz> \
    --cadd_wg_snvs_tbi </path/to/CADD/cache/whole_genome_SNVs.tsv.gz.tbi>
```

### Downloading CADD files

An helper script has been designed to help downloading CADD files.
Such files are meant to be share between multiple users, so this script is mainly meant for people administrating servers, clusters and advanced users.

```bash
nextflow run download_cache.nf --cadd_cache </path/to/CADD/cache> --cadd_version <CADD version> --genome <GENOME>
```

## Using VEP GeneSplicer plugin

To enable the use of the VEP GeneSplicer plugin:

- use the `--genesplicer` flag

Example:

```bash
nextflow run nf-core/sarek --step annotate --tools VEP --sample <file.vcf.gz> --genesplicer
```
