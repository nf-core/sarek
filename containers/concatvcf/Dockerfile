FROM debian:8.6

LABEL author="Maxime Garcia" \
description="pigz image for use in CAW" \
maintainer="maxime.garcia@scilifelab.se"

# Install libraries
RUN apt-get update && apt-get install -y --no-install-recommends \
pigz \
&& rm -rf /var/lib/apt/lists/*

# Create UPPMAX directories
RUN mkdir /pica /proj /sw
