FROM maxulysse/samtools:1.1

LABEL author="Maxime Garcia" \
description="BWA 0.7.8 and Samblaster 0.1.22 image for use in CAW" \
maintainer="maxime.garcia@scilifelab.se"

# Install libraries
RUN apt-get update && apt-get install -y --no-install-recommends \
build-essential \
ca-certificates \
curl \
git \
libncurses5-dev \
libncursesw5-dev \
python3 \
zlib1g-dev \
&& rm -rf /var/lib/apt/lists/*

# Setup ENV variables
ENV BWA_VERSION="0.7.8" \
SAMBLASTER_VERSION="v.0.1.22"

# Install BWA and SamBlaster
RUN git clone http://github.com/lh3/bwa.git bwa \
--branch $BWA_VERSION \
&& cd bwa \
&& make \
&& cp bwa /usr/local/bin/bwa \
&& cd .. \
&& git clone git://github.com/GregoryFaust/samblaster.git samblaster \
--branch $SAMBLASTER_VERSION \
&& cd samblaster \
&& make \
&& cp samblaster /usr/local/bin/samblaster \
&& rm -rf bwa samblaster
