FROM openjdk:8-slim

LABEL \
  author="Maxime Garcia" \
  description="snpEff image for use in Sarek" \
  maintainer="maxime.garcia@scilifelab.se"

# Install libraries
RUN \
  apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    libncurses5-dev \
    libncursesw5-dev \
    wget \
    zlib1g-dev \
  && rm -rf /var/lib/apt/lists/*

# Setup ENV variables
ENV \
  SNPEFF_VERSION=4_3r \
  SNPEFF_HOME=/opt/snpEff

# Install snpEff
RUN \
  wget --quiet -O snpEff_v${SNPEFF_VERSION}_core.zip \
    http://downloads.sourceforge.net/project/snpeff/snpEff_v${SNPEFF_VERSION}_core.zip \
  && unzip snpEff_v${SNPEFF_VERSION}_core.zip -d /opt/ \
  && rm snpEff_v${SNPEFF_VERSION}_core.zip
