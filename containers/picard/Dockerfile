FROM openjdk:8-slim

LABEL author="Maxime Garcia" \
description="Picard 2.0.1 image for use in CAW" \
maintainer="maxime.garcia@scilifelab.se"

# Install libraries
RUN apt-get update && apt-get install -y --no-install-recommends \
curl \
&& rm -rf /var/lib/apt/lists/*

# Setup ENV variables
ENV PICARD_BIN="picard-tools-2.0.1.zip" \
PICARD_HOME="/opt/picard-tools-2.0.1" \
PICARD_VERSION="2.0.1"

# Install PicardTools 2.0.1
RUN curl -fsSL https://github.com/broadinstitute/picard/releases/download/$PICARD_VERSION/$PICARD_BIN -o /opt/$PICARD_BIN \
&& unzip /opt/$PICARD_BIN -d /opt/ \
&& rm /opt/$PICARD_BIN

# Create UPPMAX directories
RUN mkdir /pica /proj /sw
