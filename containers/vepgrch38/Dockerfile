FROM maxulysse/vep:1.1

LABEL author="Maxime Garcia" \
description="VEP 90.1 image with GRCH38 genome for use in CAW" \
maintainer="maxime.garcia@scilifelab.se"

# Setup ENV variables
ENV GENOME="GRCh38"

# Download Genome
RUN mkdir -p $HOME/.vep \
&& cd $HOME/.vep \
&& curl -O ftp://ftp.ensembl.org/pub/release-90/variation/VEP/homo_sapiens_vep_90_$GENOME.tar.gz \
&& tar xzf homo_sapiens_vep_90_$GENOME.tar.gz \
&& rm homo_sapiens_vep_90_$GENOME.tar.gz
