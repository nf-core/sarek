# Annotation

## Tools

Within Sarek, annotation is done using snpEff, VEP, or both consecutively:
- `--tools snpEff`
  - To annotate using `snpEff`
- `--tools VEP`
  - To annotate using `VEP`
- `--tools snpEff,VEP`
  - To annotate using `snpEff` and `VEP`
- `--tools merge`
  - To annotate using `snpEff` then `VEP`

## Using containers

Sarek has already designed containers with `snpEff` and `VEP` files for `GRCh37` and `GRCh38`.
Default settings will run using these containers.

## Using cache

Both `snpEff` and `VEP` enable usage of cache.
If cache is available on the machine where Sarek is run, it is possible to run annotation using cache.
You need to specify the cache directory using `--snpEff_cache` and `--vep_cache` in the command lines or within configuration files.
The cache will only be used when `--annotation_cache` and cache directories are specified (either in command lines or in a configuration file).

Example:
```
nextflow run annotate.nf --tools snpEff --annotateVCF file.vcf.gz --genome GRCh38 --snpEff_cache /Path/To/snpEffCache --annotation_cache
nextflow run annotate.nf --tools VEP --annotateVCF file.vcf.gz --genome GRCh38 --vep_cache /Path/To/vepCache --annotation_cache
```

## Using VEP CADD plugin

To enable the use of the VEP CADD plugin:
 - Download the CADD files
 - Specify them (either on the command line, like in the example or in a configuration file)
 - use the `--cadd_cache` flag

Example:
```
nextflow run annotate.nf --tools VEP --annotateVCF file.vcf.gz --genome GRCh38 --cadd_cache \
--cadd_InDels /PathToCADD/InDels.tsv.gz \
--cadd_InDels_tbi /PathToCADD/InDels.tsv.gz.tbi \
--cadd_WG_SNVs /PathToCADD/whole_genome_SNVs.tsv.gz \
--cadd_WG_SNVs_tbi /PathToCADD/whole_genome_SNVs.tsv.gz.tbi
```

### Downloading CADD files

An helper script has been designed to help downloading CADD files.
Such files are meant to be share between multiple users, so this script is mainly meant for people administrating servers, clusters and advanced users.
```
nextflow run build.nf --cadd_cache /Path/To/CADDcache --genome <GENOME>
```

## Using VEP GeneSplicer plugin

To enable the use of the VEP GeneSplicer plugin:
 - use the `--genesplicer` flag

Example:
```
nextflow run annotate.nf --tools VEP --annotateVCF file.vcf.gz --genome GRCh38 --genesplicer
```
