FROM openjdk:8

LABEL \
  author="Maxime Garcia" \
  description="FastQC image for use in CAW" \
  maintainer="maxime.garcia@scilifelab.se"

# Install libraries
RUN \
  apt-get update && apt-get install -y --no-install-recommends \
    wget \
  && rm -rf /var/lib/apt/lists/*

# Setup ENV variables
ENV FASTQC_VERSION=0.11.5

# Install FastQC
RUN \
  wget --quiet -O fastqc_v${FASTQC_VERSION}.zip \
    http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v${FASTQC_VERSION}.zip \
  && unzip fastqc_v${FASTQC_VERSION}.zip -d /opt/ \
  && chmod 755 /opt/FastQC/fastqc \
  && ln -s /opt/FastQC/fastqc /usr/local/bin/fastqc \
  && rm fastqc_v${FASTQC_VERSION}.zip

# Create UPPMAX directories
RUN mkdir /pica /proj /scratch /sw
