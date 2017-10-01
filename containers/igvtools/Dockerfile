FROM openjdk:8-slim

LABEL author="Maxime Garcia" \
description="IGVtools 2.3.98 image for use in CAW" \
maintainer="maxime.garcia@scilifelab.se"

# Install libraries
RUN apt-get update && apt-get install -y --no-install-recommends \
curl \
&& rm -rf /var/lib/apt/lists/*

# Setup ENV variables
ENV \
  IGVTOOLS_HOME=/opt/IGVTools \
  IGVTOOLS_VERSION=2.3.98

# Install IGVTools
RUN curl -fsSL http://data.broadinstitute.org/igv/projects/downloads/2.3/igvtools_${IGVTOOLS_VERSION}.zip -o igvtools_${IGVTOOLS_VERSION}.zip \
&& unzip igvtools_${IGVTOOLS_VERSION}.zip \
&& rm igvtools_${IGVTOOLS_VERSION}.zip \
&& mv IGVTools $IGVTOOLS_HOME

RUN ls -l $IGVTOOLS_HOME

# Create UPPMAX directories
RUN mkdir /pica /proj /sw
