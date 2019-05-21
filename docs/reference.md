# Genomes and reference files

## AWS iGenomes
Sarek is using [AWS iGenomes](https://ewels.github.io/AWS-iGenomes/), which facilitate storing and sharing references.
Sarek currently uses `GRCh38` by default.
Both `GRCh37` and `GRCh38` are available with `--genome GRCh37` or `--genome GRCh38` respectively with any profile using the `conf/igenomes.config` file, or you can specify it with `-c conf/igenomes.config`.
Use `--genome smallGRCh37` to map against a small reference genome based on GRCh37.
Settings in `igenomes.config` can be tailored to your needs.

## build.nf

The [`build.nf`](#buildnf) script is used to build reference needed for smallGRCh37.

```
nextflow run build.nf
```
