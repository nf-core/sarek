FROM nfcore/base:latest

LABEL \
  author="Maxime Garcia" \
  description="snpEff image for use in nf-core/sarek" \
  maintainer="maxime.garcia@scilifelab.se"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/sarek-snpeff-2.5dev/bin:$PATH

# Setup default ARG variables
ARG CACHE_VERSION=86
ARG GENOME=GRCh38

# Download Genome
RUN snpEff download -v ${GENOME}.{CACHE_VERSION}
