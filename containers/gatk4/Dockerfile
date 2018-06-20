FROM openjdk:8-slim

LABEL \
  author="Maxime Garcia" \
  description="GATK image for use in Sarek" \
  maintainer="maxime.garcia@scilifelab.se"

# Install libraries
RUN \
  apt-get update && apt-get install -y --no-install-recommends \
  wget \
  && rm -rf /var/lib/apt/lists/*

# Setup ENV variables
ENV \
  GATK_VER="4.0.5.1" \
  GATK_HOME="/opt/gatk"

# Install GATK
RUN \
  mkdir $GATK_HOME \
  && cd $GATK_HOME \
  && wget --quiet -O $GATK_HOME/gatk-${GATK_VER}.zip \
  https://github.com/broadinstitute/gatk/releases/download/${GATK_VER}/gatk-${GATK_VER}.zip \
  && unzip gatk-${GATK_VER}.zip \
  && mv gatk-${GATK_VER}/gatk-package-${GATK_VER}-local.jar gatk4.jar \
  && rm gatk-${GATK_VER}.zip
