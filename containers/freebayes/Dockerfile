FROM debian:8.6

LABEL author="Maxime Garcia" \
description="FreeBayes 1.1.0 image for use in CAW" \
maintainer="maxime.garcia@scilifelab.se"

# Install libraries
RUN apt-get update && apt-get install -y --no-install-recommends \
build-essential \
ca-certificates \
cmake \
git \
wget \
zlib1g-dev \
&& rm -rf /var/lib/apt/lists/*

# Setup ENV variables
ENV FREEBAYES_BIN="v1.1.0.zip" \
FREEBAYES_VERSION="v1.1.0"

# Install BCFTools
RUN git clone --recursive git://github.com/ekg/freebayes.git freebayes \
--branch $FREEBAYES_VERSION \
&& cd freebayes \
&& make \
&& make install \
&& cd .. \
&& rm -rf freebayes

# Create UPPMAX directories
RUN mkdir /pica /proj /sw
