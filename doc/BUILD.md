# Building

Use the Nextflow script to build and/or push containers from Docker and/or Singularity.

All the containers have built in UPPMAX directories, so there is no need to add them for use on UPPMAX clusters.
- See the [Singularity UPPMAX guide](https://www.uppmax.uu.se/support-sv/user-guides/singularity-user-guide/)

## Usage

```bash
nextflow run . [--docker] [--singularity] [--containerPath <path>] [--push] [--containers <container1[,container2..]>] [--repository <repository>] [--tag tag]
```

- `--containers`: Choose which containers to build. Default: `all`. Possible values (to separate by commas):
  - `all` -  Build all available containers.
  - `fastqc`
  - `freebayes`
  - `gatk`
  - `igvtools`
  - `multiqc`
  - `mutect1`
  - `picard`
  - `qualimap`
  - `r-base`
  - `runallelecount`
  - `sarek`
  - `snpeff` this container serves as a base for `snpeffgrch37` and `snpeffgrch38`
  - `snpeffgrch37`
  - `snpeffgrch38`
  - `vepgrch37`
  - `vepgrch38`

- `--docker`: Build containers using `Docker`
- `--push`: Push containers to `DockerHub`
- `--repository`: Build containers under given repository. Default: `maxulysse`
- `--singularity`: Build containers using `Singularity`.
- `--containerPath`: Select where to download containers. Default: `$PWD`
- `--tag`: Build containers using given tag. Default is version number.

## Example

```bash
nextflow run . --docker --singularity --push --containers multiqc,fastqc
```

## For lazy users
We provide script to build/push or pull all containers
```bash
./scripts/do_all.sh        # Build all docker containers
./scripts/do_all.sh --push # Build and push all Docker containers into DockerHub
./scripts/do_all.sh --pull # Pull all containers from DockerHub into Singularity
```

---
[![](images/SciLifeLab_logo.png "SciLifeLab")][scilifelab-link]
[![](images/NGI_logo.png "NGI")][ngi-link]
[![](images/NBIS_logo.png "NBIS")][nbis-link]

[nbis-link]: https://www.nbis.se/
[ngi-link]: https://ngisweden.scilifelab.se/
[scilifelab-link]: https://www.scilifelab.se/
