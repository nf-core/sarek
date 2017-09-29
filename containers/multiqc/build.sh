#!/bin/sh
set -eu pipefail

# Local ENV variables
MULTIQC_VERSION=1.1

# Install libraries
apt-get update && apt-get install -y --no-install-recommends \
g++ \
git \
wget

# Install pip
wget --quiet -O /opt/get-pip.py https://bootstrap.pypa.io/get-pip.py
python /opt/get-pip.py

# Install tools
pip install networkx==1.11 # fix issue 592
pip install git+git://github.com/ewels/MultiQC.git@v${MULTIQC_VERSION}

# Clean up install
cd /
apt-get remove -y \
g++ \
git \
wget
rm -rf /build /var/lib/apt/lists/* /opt/get-pip.py

# Create UPPMAX directories
mkdir /pica /proj /sw
