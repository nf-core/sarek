# Tests

One script is available for testing purpose:
- [`scripts/test.sh`](../scripts/test.sh)

Two optional arguments are supported:
- `-p` || `--profile`:
  - `dockerTest` test using Docker containers
  - `singularityTest` (default) test using Singularity containers
- `-t` || `--test`::
 - `MAPPING`: will try preprocessing
 - `REALIGN`: will try realignment
 - `RECALIBRATE`: will try recalibration
 - `ANNOTATE`: will try variant calling and annotation
 - `ALL`: will try all the previous tests (default)

## Usage

```bash
./scripts/test.sh                # will try all tests using Singularity
./scripts/test.sh -p dockerTest  # will try all tests using Docker
./scripts/test.sh -t MAPPING     # will try MAPPING tests using Singularity
```

--------------------------------------------------------------------------------

[![](images/SciLifeLab_logo.png "SciLifeLab")][scilifelab-link]
[![](images/NGI_logo.png "NGI")][ngi-link]
[![](images/NBIS_logo.png "NBIS")][nbis-link]

[nbis-link]: https://www.nbis.se/
[ngi-link]: https://ngisweden.scilifelab.se/
[scilifelab-link]: https://www.scilifelab.se/
