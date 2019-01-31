# Configuration and Profiles

For more informations on how to use configuration files, have a look at the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html)

For more informations about profiles, have a look at the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html#config-profiles)

We provides several configuration files and profiles for Sarek.
The standard ones are designed to work on a Swedish UPPMAX cluster, but can be modified and tailored to your own need.


## Configuration files

Every configuration file can be modified for your own use.
If you want you can specify the use of a config file using `-c <config file>`

### [`aws-batch.config`](https://github.com/SciLifeLab/Sarek/blob/master/conf/aws-batch.config)

Designed for usage with AWS batch.

### [`base.config`](https://github.com/SciLifeLab/Sarek/blob/master/conf/base.config)

Define default parameters, is included into every profiles.

### [`binac.config`](https://github.com/SciLifeLab/Sarek/blob/master/conf/binac.config)

Define usage limits and Singularity for BINAC cluster in Tuebingen.

### [`cfc.config`](https://github.com/SciLifeLab/Sarek/blob/master/conf/cfc.config)

Designed for usage with Singularity on CFC at QBic.

### [`conda.config`](https://github.com/SciLifeLab/Sarek/blob/master/conf/conda.config)

> /!\\ Under development.

Define conda environement.

### [`containers.config`](https://github.com/SciLifeLab/Sarek/blob/master/conf/containers.config)

Define Containers for all process.
Images will be pulled automatically.
Use in your own profile if needed.

### [`docker.config`](https://github.com/SciLifeLab/Sarek/blob/master/conf/docker.config)

Specify Docker options.
To be used with [`containers.config`](#containersconfig)

### [`genomes.config`](https://github.com/SciLifeLab/Sarek/blob/master/conf/genomes.config)

Contain path to all references.
Modify it if you want to change genome version, or the path to your references files.

### [`igenomes.config`](https://github.com/SciLifeLab/Sarek/blob/master/conf/igenomes.config)

Contain path to all AWS iGenomes references.
Modify it if you want to change genome version, or the path to your references files.

### [`munin.config`](https://github.com/SciLifeLab/Sarek/blob/master/conf/munin.config)

Define usage limits and Singularity for munin server at BTB.

### [`resources.config`](https://github.com/SciLifeLab/Sarek/blob/master/conf/resources.config)

Define Generalized resource configuration for clusters.

### [`singularity-path.config`](https://github.com/SciLifeLab/Sarek/blob/master/conf/singularity-path.config)

Define path to Singularity Containers for all process.
To be used when downloading Singularity Containers, like on a secure UPPMAX cluster.
Images will not be pulled automatically.
You need to set them up before.

### [`singularity.config`](https://github.com/SciLifeLab/Sarek/blob/master/conf/singularity.config)

Specify Singularity options.
To be used with [`containers.config`](#containersconfig)

### [`travis.config`](https://github.com/SciLifeLab/Sarek/blob/master/conf/travis.config)

To be used for Travis (2 cpus) or on small computer for testing purpose

### [`uppmax-localhost.config`](https://github.com/SciLifeLab/Sarek/blob/master/conf/uppmax-localhost.config)

Local configuration for a UPPMAX cluster
To be run on a single node
### [`uppmax-slurm.config`](https://github.com/SciLifeLab/Sarek/blob/master/conf/uppmax-slurm.config)

Slurm configuration for a UPPMAX cluster
Will run the workflow on `/scratch` using the Nextflow [`scratch`](https://www.nextflow.io/docs/latest/process.html#scratch) directive

## Profiles
A profile is a convenient way of specifying which set of configuration files to use.
The default profile is `standard`, but Sarek has multiple predefined profiles which are listed below that can be specified by specifying `-profile <profile>`:

```bash
nextflow run SciLifeLab/Sarek --sample mysample.tsv -profile myprofile
```

### `awsbatch`

This is the profile for use with AWS Batch.

### `binac`

This is the profile for use on the german BinAC cluster.

### `btb`

This is the profile for use on the BTB server munin.

### `cfc`

This is the profile for use on the CFC cluster in Tuebingen.

### `conda`

> /!\\ Under development.

This is the profile for conda testing on a small machine, or on Travis CI.
Conda environement will be built automatically.

### `docker`

This is the profile for docker testing on a small machine, or on Travis CI.
Docker images will be pulled automatically.

### `singularity`

This is the profile for Singularity testing on a small machine, or on Travis CI.
Singularity images will be pulled automatically.

### `singularityPath`

This is the profile for Singularity testing on a small machine.
Singularity images needs to be set up.

### `slurm`

This is another profile for use on a UPPMAX cluster using the job scheduler slurm with Singularity.
Will run the workflow on `/scratch`.
Singularity images are already set up.

### `standard`

This is the default profile for use on a localhost on a UPPMAX cluster with Singularity.
Singularity images are already set up.

## Customisation
The recommended way to use custom settings is to supply Sarek with an additional configuration file. You can use the files in the [`conf/`](https://github.com/SciLifeLab/Sarek/tree/master/conf) directory as an inspiration to make this new `.config` file and specify it using the `-c` flag:

```bash
nextflow run SciLifeLab/Sarek --sample mysample.tsv -c conf/personal.config
```

Any configuration field specified in this file has precedence over the predefined configurations but any field left out from the file will be set by the normal configuration files included in the specified (or `standard`) profile.

Furthermore, to find out which configuration files take action for the different profiles, the profiles are defined in the file  [`nextflow.config`](https://github.com/SciLifeLab/Sarek/blob/master/nextflow.config).
