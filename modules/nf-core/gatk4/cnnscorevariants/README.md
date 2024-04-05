# Conda is not supported at the moment

The [bioconda](https://bioconda.github.io/recipes/gatk4/README.html) recipe is not fully working as expected, cf [github issue](https://github.com/broadinstitute/gatk/issues/7811)

Hence, we are using the docker container provided by the authors of the tool:

- [broadinstitute/gatk](https://hub.docker.com/r/broadinstitute/gatk)

This image is mirrored on the [nf-core quay.io](https://quay.io/repository/nf-core/gatk) for convenience.
