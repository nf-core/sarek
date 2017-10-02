FROM debian:8.9

LABEL \
  author="Maxime Garcia" \
  description="Image with tools used in CAW" \
  maintainer="maxime.garcia@scilifelab.se"

# Export PATH
ENV \
  MANTA_INSTALL_PATH=/opt/manta \
  PATH=/opt/bcftools/bin:/opt/htslib/bin:/opt/samtools/bin:$PATH \
  STRELKA_INSTALL_PATH=/opt/strelka

# Build container
ADD \
  build.sh /usr/bin/
RUN \
  build.sh
