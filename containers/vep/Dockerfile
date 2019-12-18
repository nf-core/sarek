FROM nfcore/base:1.7

LABEL \
  author="Maxime Garcia" \
  description="VEP image for use in nf-core/sarek" \
  maintainer="maxime.garcia@scilifelab.se"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/sarek-vep-2.5.2/bin:$PATH

# Setup default ARG variables
ARG GENOME=GRCh38
ARG SPECIES=homo_sapiens
ARG VEP_VERSION=95

# Download Genome
RUN vep_install \
  -a c \
  -c .vep \
  -s ${SPECIES} \
  -v ${VEP_VERSION} \
  -y ${GENOME} \
  --CACHE_VERSION ${VEP_VERSION} \
  --CONVERT \
  --NO_BIOPERL --NO_HTSLIB --NO_TEST --NO_UPDATE
