FROM maxulysse/samtools:1.1

LABEL author="Maxime Garcia" \
description="Manta 1.1.1 and SAMTools 1.4 image for use in CAW" \
maintainer="maxime.garcia@scilifelab.se"

# Install libraries
RUN apt-get update && apt-get install -y --no-install-recommends \
build-essential \
ca-certificates \
curl \
libncurses5-dev \
libncursesw5-dev \
python2.4 \
zlib1g-dev \
&& rm -rf /var/lib/apt/lists/*

# Setup ENV variables
ENV MANTA_BIN="manta-1.1.1.centos5_x86_64.tar.bz2" \
MANTA_INSTALL_PATH="/opt/manta" \
MANTA_VERSION="1.1.1"

# Install Manta
RUN curl -fsSL https://github.com/Illumina/manta/releases/download/v$MANTA_VERSION/$MANTA_BIN -o /opt/$MANTA_BIN \
&& tar xvjf /opt/$MANTA_BIN -C /opt/ \
&& mv /opt/manta-$MANTA_VERSION.centos5_x86_64 $MANTA_INSTALL_PATH \
&& rm /opt/$MANTA_BIN
