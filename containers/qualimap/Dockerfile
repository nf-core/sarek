FROM openjdk:8-slim

LABEL \
  author="Maxime Garcia" \
  description="Qualimap image for use in CAW" \
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

# Install BamQC + fix bugs
# See https://github.com/docker-library/tomcat/issues/80
RUN \
  rm /usr/lib/jvm/java-8-openjdk-amd64/jre/lib/accessibility.properties \
    /etc/java-8-openjdk/accessibility.properties \
  && wget --quiet -O qualimap_v${QUALIMAP_VERSION}.zip \
    https://bitbucket.org/kokonech/qualimap/downloads/qualimap_v${QUALIMAP_VERSION}.zip \
  && unzip qualimap_v${QUALIMAP_VERSION}.zip -d /opt/ \
  && rm qualimap_v${QUALIMAP_VERSION}.zip \
  && mv /opt/qualimap_v${QUALIMAP_VERSION} /opt/qualimap

# Create UPPMAX directories
RUN mkdir /pica /proj /scratch /sw
