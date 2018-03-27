FROM debian:8.9

LABEL \
  author="Maxime Garcia" \
  description="alleleCount 2.2.0 used in Sarek" \
  maintainer="maxime.garcia@scilifelab.se"

# Export PATH
ENV \
  ALLELECOUNT_VERSION=2.2.0 \
  PATH=/opt/bin:$PATH \
  PERL5LIB=$PERL5LIB:/opt/lib/perl5

# Install libraries
RUN \
  apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    libncursesw5-dev \
    wget \
    zlib1g-dev \
  && rm -rf /var/lib/apt/lists/*

# Install alleleCount
RUN \
  wget --quiet -O /opt/v${ALLELECOUNT_VERSION}.tar.gz \
    https://github.com/cancerit/alleleCount/archive/v${ALLELECOUNT_VERSION}.tar.gz \
  && tar xvzf /opt/v${ALLELECOUNT_VERSION}.tar.gz -C /opt/ \
  && cd /opt/alleleCount-${ALLELECOUNT_VERSION} \
  && ./setup.sh /opt/ \
  && rm /opt/v${ALLELECOUNT_VERSION}.tar.gz
