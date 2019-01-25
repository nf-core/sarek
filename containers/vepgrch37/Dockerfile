FROM nfcore/base:latest

LABEL \
  author="Maxime Garcia" \
  description="VEP image with GRCh37 version 95 genome for use in Sarek" \
  maintainer="maxime.garcia@scilifelab.se"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/sarek-snpeff-2.2.2/bin:$PATH

# Setup ENV variables
ENV \
  GENOME=GRCh37 \
  VEP_VERSION=95

# Download Genome
RUN \
  mkdir -p $HOME/.vep \
  && cd $HOME/.vep \
  && wget --quiet -O homo_sapiens_vep_${VEP_VERSION}_${GENOME}.tar.gz \
    ftp://ftp.ensembl.org/pub/release-${VEP_VERSION}/variation/VEP/homo_sapiens_vep_${VEP_VERSION}_${GENOME}.tar.gz \
  && tar xzf homo_sapiens_vep_${VEP_VERSION}_${GENOME}.tar.gz \
  && rm homo_sapiens_vep_${VEP_VERSION}_${GENOME}.tar.gz
