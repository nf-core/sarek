# RELEASE

> This document is for helping Sarek core developers and anyone joining the team to prepare a new release

## [CHECKLIST](https://github.com/SciLifeLab/Sarek/blob/master/.github/RELEASE_CHECKLIST.md)

This checklist is for our own reference, to help us prepare a new release.
Just follow it and be sure to check every item on the list.

## [Helper script](https://github.com/SciLifeLab/Sarek/blob/master/scripts/do_release.sh)

This script will update the version number in the following files:

-   [CHANGELOG.md](https://github.com/SciLifeLab/Sarek/blob/master/CHANGELOG.md)
    -   Will change Unreleased to correct version number and add codename and date
-   [Dockerfile](https://github.com/SciLifeLab/Sarek/blob/master/Dockerfile)
    -   Will update to correct version number
-   [Singularity](https://github.com/SciLifeLab/Sarek/blob/master/Singularity)
    -   Will update to correct version number
-   [conf/base.config](https://github.com/SciLifeLab/Sarek/blob/master/conf/base.config)
    -   Will update to correct version number

### Usage

### Usage

```bash
./scripts/do_release.sh -r "<RELEASE>" -c "<CODENAME>"
```

-   `-r|--release` specify the new version number
-   `-c|--codename` specify the codename

### Example

```bash
./scripts/do_release.sh -r "2.2.0" -c "Skårki"
```
