FROM openjdk:8-slim

LABEL author="Maxime Garcia" \
description="IGVtools 2.3.98 image for use in CAW" \
maintainer="maxime.garcia@scilifelab.se"

# Install libraries
RUN apt-get update && apt-get install -y --no-install-recommends \
curl \
&& rm -rf /var/lib/apt/lists/*

# Setup ENV variables
ENV IGVTOOLS_BIN="igvtools_2.3.98.zip" \
IGVTOOLS_HOME="/opt/IGVTools"

# Install IGVTools
RUN curl -fsSL http://data.broadinstitute.org/igv/projects/downloads/$IGVTOOLS_BIN -o /opt/$IGVTOOLS_BIN \
&& cd /opt \
&& unzip $IGVTOOLS_BIN \
&& rm $IGVTOOLS_BIN

# Create UPPMAX directories
RUN mkdir /pica /proj /sw
