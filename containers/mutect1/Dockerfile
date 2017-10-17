FROM openjdk:7-slim

LABEL \
  author="Maxime Garcia" \
  description="MuTect 1.1.5 image for use in CAW" \
  maintainer="maxime.garcia@scilifelab.se"

# Install libraries
RUN \
  apt-get update && apt-get install -y --no-install-recommends \
    wget \
  && rm -rf /var/lib/apt/lists/*

# Setup ENV variables
ENV \
  MUTECT_HOME=/opt/mutect-1.1.5 \
  MUTECT_VERSION=1.1.5

# Install MuTect1
RUN \
  wget --quiet -O muTect-${MUTECT_VERSION}-bin.zip \
    https://github.com/broadinstitute/mutect/releases/download/${MUTECT_VERSION}/muTect-${MUTECT_VERSION}-bin.zip \
  && unzip muTect-${MUTECT_VERSION}-bin.zip -d ${MUTECT_HOME} \
  && rm muTect-${MUTECT_VERSION}-bin.zip \
&& mv ${MUTECT_HOME}/muTect-${MUTECT_VERSION}.jar ${MUTECT_HOME}/muTect.jar
# Create UPPMAX directories
RUN mkdir /pica /proj /scratch /sw
