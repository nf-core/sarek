FROM openjdk:8

LABEL \
  author="Maxime Garcia" \
  description="Qualimap image for use in Sarek" \
  maintainer="maxime.garcia@scilifelab.se"

# Install libraries
RUN \
  apt-get update && apt-get install -y --no-install-recommends \
    wget \
  && rm -rf /var/lib/apt/lists/*

# Setup ENV variables
ENV \
  PATH=$PATH:/opt/qualimap \
  QUALIMAP_VERSION=2.2.1

# Install BamQC
RUN \
  wget --quiet -O qualimap_v${QUALIMAP_VERSION}.zip \
    https://bitbucket.org/kokonech/qualimap/downloads/qualimap_v${QUALIMAP_VERSION}.zip \
  && unzip qualimap_v${QUALIMAP_VERSION}.zip -d /opt/ \
  && rm qualimap_v${QUALIMAP_VERSION}.zip \
  && mv /opt/qualimap_v${QUALIMAP_VERSION} /opt/qualimap

# Create UPPMAX directories
RUN mkdir /pica /proj /scratch /sw
