FROM nfcore/base:latest

LABEL \
  author="Maxime Garcia" \
  description="VEP image with GRCh38 version 95 genome for use in Sarek" \
  maintainer="maxime.garcia@scilifelab.se"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/sarek-vep-2.2.2/bin:$PATH

# Setup ENV variables
ENV \
  GENOME=GRCh38 \
  VEP_VERSION=95

# Download Genome
RUN mkdir -p .vep && cd .vep \
  && wget --quiet -O homo_sapiens_vep_${VEP_VERSION}_${GENOME}.tar.gz \
    ftp://ftp.ensembl.org/pub/release-${VEP_VERSION}/variation/VEP/homo_sapiens_vep_${VEP_VERSION}_${GENOME}.tar.gz \
  && tar xzf homo_sapiens_vep_${VEP_VERSION}_${GENOME}.tar.gz \
  && rm -rf homo_sapiens_vep_${VEP_VERSION}_${GENOME}.tar.gz
