FROM debian:8.6

LABEL author="Maxime Garcia" \
description="Strelka 2.8.2 image for use in CAW" \
maintainer="maxime.garcia@scilifelab.se"

# Install libraries
RUN apt-get update && apt-get install -y --no-install-recommends \
bzip2 \
ca-certificates \
curl \
python \
&& rm -rf /var/lib/apt/lists/*

# Setup ENV variables
ENV STRELKA_BIN="strelka-2.8.2.centos5_x86_64.tar.bz2" \
STRELKA_FOLDER="strelka-2.8.2.centos5_x86_64" \
STRELKA_INSTALL_DIR="/opt/strelka" \
STRELKA_VERSION="2.8.2"

# Install Strelka
RUN curl -fsSL https://github.com/Illumina/strelka/releases/download/v$STRELKA_VERSION/$STRELKA_BIN -o $STRELKA_BIN \
&& tar -xjf $STRELKA_BIN \
&& mv $STRELKA_FOLDER $STRELKA_INSTALL_DIR \
&& rm -rf $STRELKA_BIN*

# Create UPPMAX directories
RUN mkdir /pica /proj /sw
