# Tests

One script is available for testing purpose:
- [`scripts/test.sh`](../scripts/test.sh)

Four optional arguments are supported:
- `-g` || `--genome`:
  Choose the genome reference version (overwrite configuration files and profiles)
- `-p` || `--profile`:
  Choose which profile to test. These options should work on a personnal computer:
  - `docker` test using Docker containers
  - `singularity` (default) test using Singularity containers
- `-s` || `--sample`:
  Use to change the test sample (default=`data/tsv/tiny.tsv`)
- `-t` || `--test`:
 - `DIR`: test `mapping` with an input directory, all other tests use a TSV file
 - `STEP`: test `mapping`, `realign` and `recalibrate`
 - `GERMLINE`: test `mapping` and Variant Calling with `HaplotypeCaller`
 - `TOOLS`: test `mapping` and Variant Calling with `FreeBayes`, `HaplotypeCaller`, `MuTect1`, `MuTect2`, `Strelka`
 - `MANTA`: test `mapping` and Variant Calling with `Manta`
 - `ANNOTATESNPEFF`: test annotation using `snpEFF`
 - `ANNOTATEVEP`: test annotation using `VEP`
 - `BUILDCONTAINERS`: test building all containers except `snpeffgrch37`, `snpeffgrch38`, `vepgrch37` and `vepgrch38`
 - `ALL`: test all the previous tests (default)

## Usage

```bash
# Will try all tests using Singularity
./scripts/test.sh
# Will try all tests using Docker
./scripts/test.sh -p docker
# Will try `STEP` tests using Singularity
./scripts/test.sh -t `STEP`
# Will try `STEP` tests using Singularity with GRCh37 genome
./scripts/test.sh -t `STEP` -g GRCh37
# Will try all tests using Singularity on manta test data
./scripts/test.sh -s data/tsv/tiny-manta.tsv
```

--------------------------------------------------------------------------------

[![](images/SciLifeLab_logo.png "SciLifeLab")][scilifelab-link]
[![](images/NGI_logo.png "NGI")][ngi-link]
[![](images/NBIS_logo.png "NBIS")][nbis-link]

[nbis-link]: https://www.nbis.se/
[ngi-link]: https://ngisweden.scilifelab.se/
[scilifelab-link]: https://www.scilifelab.se/
