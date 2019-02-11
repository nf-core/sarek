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

## Downloading cache

An helper script has been designed to help downloading snpEff and VEP cache.
Cache is meant to be share between multiple users, so this script is mainly meant for people administrating servers, clusters and advanced users.
```
nextflow run buildReferences.nf --snpEff_cache /Path/To/snpEffCache --vep_cache /Path/To/vepCache --genome <GENOME>
```
