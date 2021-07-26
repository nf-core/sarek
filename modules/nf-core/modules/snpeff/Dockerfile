<<<<<<< HEAD:containers/snpeff/Dockerfile
FROM nfcore/base:1.10.2
=======
FROM nfcore/base:1.14
>>>>>>> 2fae194272556a8b9123813203a74ee204f7f36a:modules/nf-core/modules/snpeff/Dockerfile
LABEL \
    author="Maxime Garcia" \
    description="snpEff image for nf-core pipelines" \
    maintainer="maxime.garcia@scilifelab.se"

# Install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
<<<<<<< HEAD:containers/snpeff/Dockerfile
ENV PATH /opt/conda/envs/nf-core-sarek-snpeff-3.0dev/bin:$PATH
=======
ENV PATH /opt/conda/envs/nf-core-snpeff-5.0/bin:$PATH
>>>>>>> 2fae194272556a8b9123813203a74ee204f7f36a:modules/nf-core/modules/snpeff/Dockerfile

# Setup default ARG variables
ARG GENOME=GRCh38
ARG SNPEFF_CACHE_VERSION=99

# Download Genome
RUN snpEff download -v ${GENOME}.${SNPEFF_CACHE_VERSION}

# Dump the details of the installed packages to a file for posterity
<<<<<<< HEAD:containers/snpeff/Dockerfile
RUN conda env export --name nf-core-sarek-snpeff-3.0dev > nf-core-sarek-snpeff-3.0dev.yml
=======
RUN conda env export --name nf-core-snpeff-5.0 > nf-core-snpeff-5.0.yml
>>>>>>> 2fae194272556a8b9123813203a74ee204f7f36a:modules/nf-core/modules/snpeff/Dockerfile
