FROM openjdk:8-slim

LABEL author="Maxime Garcia" \
description="snpEff 4.3m image for use in CAW" \
maintainer="maxime.garcia@scilifelab.se"

# Install libraries
RUN apt-get update && apt-get install -y --no-install-recommends \
curl \
&& rm -rf /var/lib/apt/lists/*

# Install libraries
RUN apt-get update && apt-get install -y --no-install-recommends \
build-essential \
ca-certificates \
libncurses5-dev \
libncursesw5-dev \
wget \
zlib1g-dev \
&& rm -rf /var/lib/apt/lists/*

# Setup ENV variables
ENV SNPEFF_BIN="snpEff_v4_3m_core.zip" \
SNPEFF_HOME="/opt/snpEff"

# Install snpEff
RUN wget http://downloads.sourceforge.net/project/snpeff/$SNPEFF_BIN -O /opt/$SNPEFF_BIN \
&& unzip /opt/$SNPEFF_BIN -d /opt/ \
&& rm /opt/$SNPEFF_BIN

# Create UPPMAX directories
RUN mkdir /pica /proj /sw
