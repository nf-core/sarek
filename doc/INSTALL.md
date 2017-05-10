# Installation

To use this pipeline, you need to have a working version of Nextflow installed, and References files.
You can use a small reference genome as testing.
Nextflow can also be use conjointly with Docker to facilitate the use of other tools.

- See the [Install Nextflow documentation](https://github.com/SciLifeLab/NGI-NextflowDocs/blob/master/docs/INSTALL.md)
- See the [Reference files documentation](REFERENCES.md)

This small tutorial will explain to you how to run CAW on a small sample test data on a Swedish UPPMAX cluster. Some variables are specific, but it can be easily modified to suit any clusters.

## Make a test directory

```bash
mkdir test_CAW
cd test_CAW
```
## Build the smallGRCh37 Reference

```bash
nextflow run SciLifeLab/CAW/buildReferences.nf --download --genome smallGRCh37
```

```bash
nextflow run SciLifeLab/CAW/buildReferences.nf --refDir <pathToSmallRefRepo> --genome smallGRCh37
```
- See the [genomes files documentation](GENOMES.md)

## Test the workflow on a test tiny set

This workflow itself needs no installation. Nextflow will automatically fetch it from GitHub when launched if `SciLifeLab/CAW` is specified as the workflow name.

```bash
nextflow run SciLifeLab/CAW --test --genome smallGRCh37
```

If you're using a Swedish UPPMAX cluster, don't forget to provide your project ID.

```bash
nextflow run SciLifeLab/CAW --test --genome smallGRCh37 --project <UPPMAX Project ID>
```

# Other possibilities for advanced users

## Clone the repository and test the workflow on a test tiny set

You can download the repository yourself from GitHub and run them directly:

```bash
git clone https://github.com/SciLifeLab/CAW test_CAW
git clone https://github.com/szilvajuhos/smallRef smallGRCh37
cd test_CAW
nextflow run SciLifeLab/CAW/buildReferences.nf --refDir ../smallGRCh37 --genome smallGRCh37
nextflow run main.nf --test --genome smallGRCh37
```

## Load Nextflow

If you're running on a Swedish UPPMAX cluster you can load Nextflow as an environment module:

```bash
module load bioinfo-tools Nextflow
```

Environnement variables are set up each time the module is loaded, so you might want to set them up after loading the module.

## Running tests in interactive mode on milou

You can try the test data by changing to the interactive mode on milou and run the test tiny set like:

```
$ interactive -A <UPPMAX Project ID> -p node
[ ... login messages ... ]
$ nextflow run main.nf -profile interactive --test --project <UPPMAX Project ID>
```

For more tests, see [Usage documentation](USAGE.md#test)

--------------------------------------------------------------------------------

[![](images/SciLifeLab_logo.png "SciLifeLab")][scilifelab-link] [![](images/NGI-final-small.png "NGI")][ngi-link]

[ngi-link]: https://ngisweden.scilifelab.se/
[scilifelab-link]: http://www.scilifelab.se/
