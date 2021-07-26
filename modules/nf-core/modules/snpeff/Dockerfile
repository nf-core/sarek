FROM nfcore/base:1.14
LABEL \
    author="Maxime Garcia" \
    description="snpEff image for nf-core pipelines" \
    maintainer="maxime.garcia@scilifelab.se"

# Install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-snpeff-5.0/bin:$PATH

# Setup default ARG variables
ARG GENOME=GRCh38
ARG SNPEFF_CACHE_VERSION=99

# Download Genome
RUN snpEff download -v ${GENOME}.${SNPEFF_CACHE_VERSION}

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-snpeff-5.0 > nf-core-snpeff-5.0.yml
