# Genomes and reference files

## AWS iGenomes
Sarek is using [AWS iGenomes](https://ewels.github.io/AWS-iGenomes/), which facilitate storing and sharing references.
Both `GRCh37` and `GRCh38` are available with `--genome GRCh37` or `--genome GRCh38` respectively with any profile using the `conf/igenomes.config` file, or you can specify it with `-c conf/igenomes.config`.

Sarek currently uses `GRCh38` by default.

Settings in `igenomes.config` can be tailored to your needs.

The [`build.nf`](#buildnf) script is used to build the indexes for the reference test.

Use `--genome smallGRCh37` to map against a small reference genome based on GRCh37.

## build.nf

The `build.nf` script can build the files needed for smallGRCh37.

```
nextflow run build.nf
```
