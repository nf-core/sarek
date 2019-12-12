FROM nfcore/base:1.7

LABEL \
  author="Maxime Garcia" \
  description="snpEff image for use in nf-core/sarek" \
  maintainer="maxime.garcia@scilifelab.se"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/sarek-snpeff-2.5.2/bin:$PATH

# Setup default ARG variables
ARG GENOME=GRCh38
ARG SNPEFF_CACHE_VERSION=86

# Download Genome
RUN snpEff download -v ${GENOME}.${SNPEFF_CACHE_VERSION}
