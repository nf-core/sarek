# Configuration and Profiles

For more informations on how to use configuration files, have a look at the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html)

For more informations about profiles, have a look at the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html#config-profiles)

We provides several configuration files and profiles for Sarek. The standard ones are designed to work on a Swedish UPPMAX clusters, and can be modified and tailored to your own need.

## Configuration files

Every configuration file can be modified for your own use. If you want you can specify the use of a config file using `-c <config file>`

### [`containers.config`](https://github.com/SciLifeLab/Sarek/blob/master/configuration/containers.config)

Define Containers for all process.
Images will be pulled automatically.
Use in your own profile if needed.

### [`docker.config`](https://github.com/SciLifeLab/Sarek/blob/master/configuration/docker.config)

Define Docker Containers for all process.
Images will be pulled automatically.
Use in your own profile if needed.

### [`genomes.config`](https://github.com/SciLifeLab/Sarek/blob/master/configuration/genomes.config)

Contain path to all references.
Modify it if you want to change genome version, or the path to your references files.

### [`singularity-path.config`](https://github.com/SciLifeLab/Sarek/blob/master/configuration/singularity-path.config)

Define path to Singularity Containers for all process.
To be used when downloading Singularity Containers, like on a secure UPPMAX cluster.
Images will not be pulled automatically.
You need to set them up before.

### [`singularity.config`](https://github.com/SciLifeLab/Sarek/blob/master/configuration/singularity.config)

Define Singularity Containers for all process.
Images will be pulled automatically.
Use in your own profile if needed.

### [`travis.config`](https://github.com/SciLifeLab/Sarek/blob/master/configuration/travis.config)

To be used for Travis (2 cpus) or on small computer for testing purpose

### [`uppmax-slurm.config`](https://github.com/SciLifeLab/Sarek/blob/master/configuration/uppmax-slurm.config)

Slurm configuration for a UPPMAX cluster
Will run the workflow on `/scratch` using the Nextflow [`scratch`](https://www.nextflow.io/docs/latest/process.html#scratch) directive

## profiles

Every profile can be modified for your own use. To use a profile, you'll need to specify `-profile <profile>`

### `docker`

This is the profile for docker testing on a small machine, or on Travis CI.
Docker images will be pulled automatically.

### `standard`

This is the default profile for use on a localhost on a UPPMAX cluster with Singularity.
Singularity images need to be set up.

### `slurm`

This is another profile for use on a UPPMAX cluster using the job scheduler slurm with Singularity.
Will run the workflow on `/scratch`.
Singularity images need to be set up.

### `slurmDownload`

This is another profile for use on a UPPMAX cluster using the job scheduler slurm with Singularity.
Will run the workflow on `/scratch`.
Singularity images will be pulled automatically.

### `singularity`

This is the profile for Singularity testing on a small machine, or on Travis CI.
Singularity images will be pulled automatically.

--------------------------------------------------------------------------------

[![](images/SciLifeLab_logo.png "SciLifeLab")][scilifelab-link]
[![](images/NGI_logo.png "NGI")][ngi-link]
[![](images/NBIS_logo.png "NBIS")][nbis-link]

[nbis-link]: https://www.nbis.se/
[ngi-link]: https://ngisweden.scilifelab.se/
[scilifelab-link]: https://www.scilifelab.se/
