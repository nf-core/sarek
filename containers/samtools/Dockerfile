FROM debian:8.6

LABEL author="Maxime Garcia" \
description="SAMTools 1.4 image for use in CAW" \
maintainer="maxime.garcia@scilifelab.se"

# Install libraries
RUN apt-get update && apt-get install -y --no-install-recommends \
build-essential \
ca-certificates \
curl \
libbz2-dev \
liblzma-dev \
libncurses5-dev \
libncursesw5-dev \
zlib1g-dev \
&& rm -rf /var/lib/apt/lists/*

# Setup ENV variables
ENV SAMTOOLS_BIN="samtools-1.4.tar.bz2" \
SAMTOOLS_VERSION="1.4"

# Install SAMTools
RUN curl -fsSL https://github.com/samtools/samtools/releases/download/$SAMTOOLS_VERSION/$SAMTOOLS_BIN -o /opt/$SAMTOOLS_BIN \
&& tar xvjf /opt/$SAMTOOLS_BIN -C /opt/ \
&& cd /opt/samtools-$SAMTOOLS_VERSION \
&& make \
&& make install \
&& rm /opt/$SAMTOOLS_BIN

# Create UPPMAX directories
RUN mkdir /pica /proj /sw
