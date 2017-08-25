# Tests

One script is available for testing purpose:
- (`scripts/test.sh`)[../scripts/test.sh]

Two optional positional arguments are supported:
1. PROFILE - should be `travis` or `singularityTest` (default)
To respectively test using Docker or Singularity containers.
2. TEST - to choose between:
 - `MAPPING`: will try preprocessing
 - `REALIGN`: will try realignment
 - `RECALIBRATE`: will try recalibration
 - `ANNOTATE`: will try variant calling and annotation
 - `ALL`: will try all the previous tests (default)

## Usage

```bash
./scripts/test.sh                         # will try all tests using Singularity
./scripts/test.sh travis                  # will try all tests using Docker
./scripts/test.sh singularityTest MAPPING # will try MAPPING tests using Singularity
```

--------------------------------------------------------------------------------

[![](images/SciLifeLab_logo.png "SciLifeLab")][scilifelab-link]
[![](images/NGI_logo.png "NGI")][ngi-link]
[![](images/NBIS_logo.png "NBIS")][nbis-link]

[nbis-link]: https://www.nbis.se/
[ngi-link]: https://ngisweden.scilifelab.se/
[scilifelab-link]: https://www.scilifelab.se/
