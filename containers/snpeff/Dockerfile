FROM nfcore/base:1.9
LABEL \
  author="Maxime Garcia" \
  description="snpEff image for use in nf-core/sarek" \
  maintainer="maxime.garcia@scilifelab.se"

# Install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-sarek-snpeff-2.6dev/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-sarek-snpeff-2.6dev > nf-core-sarek-snpeff-2.6dev.yml

# Setup default ARG variables
ARG GENOME=GRCh38
ARG SNPEFF_CACHE_VERSION=86

# Download Genome
RUN snpEff download -v ${GENOME}.${SNPEFF_CACHE_VERSION}
