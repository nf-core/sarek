FROM willmclaren/ensembl-vep:release_92

LABEL \
  author="Maxime Garcia" \
  description="VEP image with GRCh38 genome for use in Sarek" \
  maintainer="maxime.garcia@scilifelab.se"

# Setup ENV variables
ENV \
  GENOME=GRCh38 \
  VEP_VERSION=91

# Download Genome
RUN \
  mkdir -p $HOME/.vep \
  && cd $HOME/.vep \
  && wget --quiet -O homo_sapiens_vep_${VEP_VERSION}_${GENOME}.tar.gz \
    ftp://ftp.ensembl.org/pub/release-${VEP_VERSION}/variation/VEP/homo_sapiens_vep_${VEP_VERSION}_${GENOME}.tar.gz \
  && tar xzf homo_sapiens_vep_${VEP_VERSION}_${GENOME}.tar.gz \
  && rm homo_sapiens_vep_${VEP_VERSION}_${GENOME}.tar.gz
