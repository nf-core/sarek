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
 - `MAPPING`: will try preprocessing
 - `REALIGN`: will try realignment
 - `RECALIBRATE`: will try recalibration
 - `ANNOTATESNPEFF`: will try variant calling and annotation using snpEff
 - `ANNOTATEVEP`: will try variant calling and annotation using VEP
 - `ALL`: will try all the previous tests (default)

## Usage

```bash
# Will try all tests using Singularity
./scripts/test.sh
# Will try all tests using Docker
./scripts/test.sh -p docker
# Will try MAPPING tests using Singularity
./scripts/test.sh -t MAPPING
# Will try MAPPING tests using Singularity with GRCh37 genome
./scripts/test.sh -t MAPPING -g GRCh37
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
