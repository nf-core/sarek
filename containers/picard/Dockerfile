FROM openjdk:8-slim

LABEL \
  author="Maxime Garcia" \
  description="Picard image for use in Sarek" \
  maintainer="maxime.garcia@scilifelab.se"

# Install libraries
RUN \
  apt-get update && apt-get install -y --no-install-recommends \
    wget \
  && rm -rf /var/lib/apt/lists/*

# Setup ENV variables
ENV \
  PICARD_HOME=/opt/picard \
  PICARD_VERSION=2.0.1

# Install PicardTools
RUN \
  wget --quiet -O picard-tools-${PICARD_VERSION}.zip \
    https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard-tools-${PICARD_VERSION}.zip \
  && unzip picard-tools-${PICARD_VERSION}.zip \
  && mv picard-tools-${PICARD_VERSION} ${PICARD_HOME} \
  && rm picard-tools-${PICARD_VERSION}.zip

# Create UPPMAX directories
RUN mkdir /pica /proj /scratch /sw
