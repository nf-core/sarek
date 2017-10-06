FROM openjdk:8-slim

LABEL \
  author="Maxime Garcia" \
  description="Picard image for use in CAW" \
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
  wget --quiet -O picard.jar https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard.jar \
  && mkdir ${PICARD_HOME} \
&& mv picard.jar ${PICARD_HOME}/picard.jar

# Create UPPMAX directories
RUN mkdir /pica /proj /scratch /sw
