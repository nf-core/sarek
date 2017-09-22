FROM debian:8.6

LABEL author="Maxime Garcia" \
description="BCFTools 1.4 image for use in CAW" \
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
ENV BCFTOOLS_BIN="bcftools-1.4.tar.bz2" \
BCFTOOLS_PLUGINS="/usr/local/libexec/bcftools" \
BCFTOOLS_VERSION="1.4" \
PATH="$PATH:/opt/bcftools/bin"

# Install BCFTools
RUN curl -fsSL https://github.com/samtools/bcftools/releases/download/$BCFTOOLS_VERSION/$BCFTOOLS_BIN -o /opt/$BCFTOOLS_BIN \
&& tar xvjf /opt/$BCFTOOLS_BIN -C /opt/ \
&& cd /opt/bcftools-$BCFTOOLS_VERSION \
&& make \
&& make prefix=/opt/bcftools install \
&& rm -rf /opt/$BCFTOOLS_BIN /opt/bcftools-$BCFTOOLS_VERSION

# Create UPPMAX directories
RUN mkdir /pica /proj /sw
