# Configuration and Profiles

For more informations on how to use configuration files, have a look at the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html)

For more informations about profiles, have a look at the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html#config-profiles)

We provides several configuration files and profiles for CAW. The standard ones are designed to work on a Swedish UPPMAX clusters, and can be modified and tailored to your own need.

## Configuration files

Every configuration file can be modified for your own use. If you want you can specify the use of a config file using `-c <config file>`

### [`docker.config`](../configuration/docker.config)

Contain Docker images for all process.
Images will be pulled automatically.
Use in your own profile if needed.

### [`genomes.config`](../configuration/genomes.config)

Contain path to all references.
Modify it if you want to change genome version, or the path to your references files.

### [`singularity-path.config`](../configuration/singularity-path.config)

To be used when downloading singularity containers, like on a secure UPPMAX cluster.
Images will not be pulled automatically.
You need to set them up before.

### [`singularity.config`](../configuration/singularity.config)

Contain Singularity images for all process.
Images will be pulled automatically.
Use in your own profile if needed.

### [`travis.config`](../configuration/travis.config)

To be used for Travis (2 cpus) or on small computer for testing purpose

### [`uppmax-localhost.config`](../configuration/uppmax-localhost.config)

To be used on a typical localhost on a UPPMAX cluster (16 cpus)

### [`uppmax-slurm.config`](../configuration/uppmax-slurm.config)

Slurm configuration for a UPPMAX cluster

## profiles

Every profile can be modified for your own use. To use a profile, you'll need to specify `-profile <profile>`

### `docker`

This is the profile for docker testing on a small machine, or on Travis CI.
Images will be pulled automatically.

### `standard`

This is the default profile for use on a localhost on a UPPMAX cluster with Singularity.
Singularity images need to be set up.

### `download`

This is the default profile for use on a localhost on a UPPMAX cluster with Singularity.
Singularity will be pulled automatically.

### `slurm`

This is another profile for use on a UPPMAX cluster using the job scheduler slurm with Singularity.
Singularity images need to be set up.

### `slurm-download`

This is another profile for use on a UPPMAX cluster using the job scheduler slurm with Singularity.
Singularity will be pulled automatically.

### `singularity`

This is the profile for Singularity testing on a small machine, or on Travis CI.
Images will be pulled automatically.

--------------------------------------------------------------------------------

[![](images/SciLifeLab_logo.png "SciLifeLab")][scilifelab-link]
[![](images/NGI_logo.png "NGI")][ngi-link]
[![](images/NBIS_logo.png "NBIS")][nbis-link]

[nbis-link]: https://www.nbis.se/
[ngi-link]: https://ngisweden.scilifelab.se/
[scilifelab-link]: https://www.scilifelab.se/
