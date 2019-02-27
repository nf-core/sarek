FROM nfcore/base:latest

LABEL \
  author="Maxime Garcia" \
  description="VEP image with GRCh38 version 95 genome for use in Sarek" \
  maintainer="maxime.garcia@scilifelab.se"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/sarek-vep-2.3/bin:$PATH

# Setup ENV variables
ENV \
  GENOME=GRCh38 \
  VEP_VERSION=95

# Download Genome
RUN vep_install \
  -a c \
  -c .vep \
  -s homo_sapiens \
  -v ${VEP_VERSION} \
  -y ${GENOME} \
  --CACHE_VERSION ${VEP_VERSION} \
  --CONVERT \
  --NO_HTSLIB --NO_TEST --NO_BIOPERL --NO_UPDATE
