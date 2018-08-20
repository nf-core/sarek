# Testing Sarek

For testing purpose we provide [Sarek-data](https://github.com/SciLifeLab/Sarek-data), a repository with test data and corresponding reference files.

One simple bash script is available, which will pull the Sarek-data repository and perform all the tests:
- [`scripts/test.sh`](../scripts/test.sh)

Such tests are used in our Continuous Integration with Travis. You can perform the same tests to familiarize yourself with the workflow.

## Testing with Singularity
For testing with Docker, just replace `singularity` with `docker` in every occurence.
```bash
# Dowload Sarek and the test data
git clone --recursive https://github.com/SciLifeLab/Sarek Sarek-test
cd Sarek-test

# Build the references for the test data
nextflow run buildReferences.nf --outDir References/smallGRCh37 \
  --refDir Sarek-data/reference --genome smallGRCh37 --tag latest \
  --verbose -profile singularity

# Testing --sampleDir as input for Germline
nextflow run main.nf --sampleDir Sarek-data/testdata/manta/normal \
  --step mapping --genome smallGRCh37 --genome_base References/smallGRCh37 \
  --tag latest -profile singularity

# Testing to restart from `realign`
nextflow run main.nf --step realign \
  --genome smallGRCh37 --genome_base References/smallGRCh37 \
  --tag latest -profile singularity

# Testing to restart from `recalibrate`
nextflow run main.nf --step recalibrate \
  --genome smallGRCh37 --genome_base References/smallGRCh37 \
  --tag latest -profile singularity

# Testing germline variant calling
nextflow run germlineVC.nf --genome smallGRCh37 \
  --genome_base References/smallGRCh37 --tools HaplotypeCaller,Manta,Strelka \
  --tag latest -profile singularity

# Testing generating report
nextflow run runMultiQC.nf -profile singularity

# Removing test data before new tests
rm -rf Preprocessing Reports VariantCalling

# Testing --sample as input for Somatic
nextflow run main.nf --sample Sarek-data/testdata/tsv/tiny-manta.tsv \
  --step mapping --genome smallGRCh37 --genome_base References/smallGRCh37 \
  --tag latest -profile singularity

# Testing germline variant calling
nextflow run germlineVC.nf --genome smallGRCh37 \
  --genome_base References/smallGRCh37 --tools HaplotypeCaller,Manta,Strelka \
  --tag latest -profile singularity

# Testing somatic variant calling
nextflow run somaticVC.nf --genome smallGRCh37 \
  --genome_base References/smallGRCh37 --tools Manta,Strelka,FreeBayes,MuTect2 \
  --tag latest -profile singularity

# Testing somatic variant calling following Strelka2 Best Practices
nextflow run somaticVC.nf --genome smallGRCh37 \
  --genome_base References/smallGRCh37 --tools Manta,Strelka,FreeBayes,MuTect2 \
  --strelkaBP --tag latest -profile singularity

# Testing annotation
nextflow run annotate.nf --tools snpEFF,VEP \
  --annotateVCF VariantCalling/StrelkaBP/Strelka_9876T_vs_1234N_somatic_indels.vcf.gz \
  -profile singularity

# Testing generating report
nextflow run runMultiQC.nf -profile singularity
```

## Testing on a secure cluster
On a secure cluster as bianca, with no internet access, you will need to download and transfer Sarek and the test data first.

Follow the [installation guide for `bianca`](https://github.com/SciLifeLab/Sarek/blob/master/doc/INSTALL_BIANCA.md).

And then start the test at the `Build the references for the test data` step.

## Usage

Four optional arguments are supported:
- `-g` || `--genome`:
  Choose the genome reference version (overwrite configuration files and profiles)
- `-p` || `--profile`:
  Choose which profile to test. These options should work on a personal computer:
  - `docker` test using Docker containers
  - `singularity` (default) test using Singularity containers
- `-s` || `--sample`:
  Use to change the test sample (default=`Sarek-data/testdata/tsv/tiny.tsv`)
- `-t` || `--test`:
 - `DIR`: test `mapping` with an input directory
 - `STEP`: test `mapping`, `realign` and `recalibrate`
 - `GERMLINE`: test `mapping` and Variant Calling with `HaplotypeCaller`
 - `TOOLS`: test `mapping` and Variant Calling with `FreeBayes`, `HaplotypeCaller`, `MuTect1`, `MuTect2`, `Strelka`
 - `MANTA`: test `mapping` and Variant Calling with `Manta`
 - `ANNOTATESNPEFF`: test annotation using `snpEFF`
 - `ANNOTATEVEP`: test annotation using `VEP`
 - `BUILDCONTAINERS`: test building all containers except `snpeffgrch37`, `snpeffgrch38`, `vepgrch37` and `vepgrch38`
 - `ALL`: test all the previous tests (default)

```bash
# Will perform all tests using Singularity
./scripts/test.sh
# Will perform all tests using Docker
./scripts/test.sh -p docker
# Will perform `STEP` tests using Singularity
./scripts/test.sh -t `STEP`
# Will perform `STEP` tests using Singularity with GRCh37 genome
./scripts/test.sh -t `STEP` -g GRCh37
# Will perform all tests using Singularity on manta test data
./scripts/test.sh -s Sarek-data/testdata/tsv/tiny-manta.tsv
```

--------------------------------------------------------------------------------

[![](images/SciLifeLab_logo.png "SciLifeLab")][scilifelab-link]
[![](images/NGI_logo.png "NGI")][ngi-link]
[![](images/NBIS_logo.png "NBIS")][nbis-link]

[nbis-link]: https://www.nbis.se/
[ngi-link]: https://ngisweden.scilifelab.se/
[scilifelab-link]: https://www.scilifelab.se/
