#!/bin/sh
set -eu pipefail

# Local ENV variables
BCFTOOLS_VERSION=1.5
BWA_VERSION=0.7.16
HTSLIB_VERSION=1.5
MANTA_VERSION=1.1.1
SAMTOOLS_VERSION=1.5
STRELKA_VERSION=2.8.2

# Install libraries
apt-get update && apt-get install -y --no-install-recommends \
  build-essential \
  bzip2 \
  ca-certificates \
  g++ \
  gcc \
  git \
  libbz2-dev \
  liblzma-dev \
  libncurses5-dev \
  libncursesw5-dev \
  make \
  python \
  python3 \
  unzip \
  wget \
  zlib1g-dev

# Install tools
mkdir /build

# Install BWA
cd /build
git clone http://github.com/lh3/bwa.git bwa \
  --branch v${BWA_VERSION}
cd bwa
make
cp bwa /usr/local/bin/bwa

# Install HTSlib
cd /build
wget --quiet -O htslib-${HTSLIB_VERSION}.tar.bz2 \
  https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2
tar xfj htslib-${HTSLIB_VERSION}.tar.bz2
rm htslib-${HTSLIB_VERSION}.tar.bz2
cd htslib-${HTSLIB_VERSION}
./configure --prefix=/opt/htslib
make && make install

# Install BCFtools
cd /build
wget --quiet -O bcftools-${BCFTOOLS_VERSION}.tar.bz2 \
  https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2
tar xfj bcftools-${BCFTOOLS_VERSION}.tar.bz2
rm bcftools-${BCFTOOLS_VERSION}.tar.bz2
cd bcftools-${BCFTOOLS_VERSION}
./configure --prefix=/opt/bcftools
make && make install

# Install Samtools
cd /build
wget --quiet -O samtools-${SAMTOOLS_VERSION}.tar.bz2 \
  https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2
tar xfj samtools-${SAMTOOLS_VERSION}.tar.bz2
rm samtools-${SAMTOOLS_VERSION}.tar.bz2
cd samtools-${SAMTOOLS_VERSION}
./configure --prefix=/opt/samtools
make && make install

# Install Manta
cd /build
wget --quiet -O manta-${MANTA_VERSION}.centos5_x86_64.tar.bz2 \
  https://github.com/Illumina/manta/releases/download/v${MANTA_VERSION}/manta-${MANTA_VERSION}.centos5_x86_64.tar.bz2
tar xvjf manta-${MANTA_VERSION}.centos5_x86_64.tar.bz2
mv manta-${MANTA_VERSION}.centos5_x86_64 $MANTA_INSTALL_PATH
rm manta-${MANTA_VERSION}.centos5_x86_64.tar.bz2

# Install Strelka
cd /build
wget --quiet -O strelka-${STRELKA_VERSION}.centos5_x86_64.tar.bz2 \
  https://github.com/Illumina/strelka/releases/download/v${STRELKA_VERSION}/strelka-${STRELKA_VERSION}.centos5_x86_64.tar.bz2
tar xvjf strelka-${STRELKA_VERSION}.centos5_x86_64.tar.bz2
mv strelka-${STRELKA_VERSION}.centos5_x86_64 $STRELKA_INSTALL_PATH
rm strelka-${STRELKA_VERSION}.centos5_x86_64.tar.bz2

# Clean up install
cd /
apt-get remove -y \
  build-essential \
  ca-certificates \
  gcc \
  git \
  libbz2-dev \
  liblzma-dev \
  libncurses5-dev \
  libncursesw5-dev \
  unzip \
  wget \
  zlib1g-dev
apt-get clean
rm -rf /build /var/lib/apt/lists/* /opt/get-pip.py
