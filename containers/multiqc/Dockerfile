FROM openjdk:8-slim

LABEL author="Maxime Garcia" \
description="MultiQC 1.0 image for use in CAW" \
maintainer="maxime.garcia@scilifelab.se"

# Install libraries
RUN apt-get update && apt-get install -y --no-install-recommends \
curl \
git \
gcc \
python-dev \
&& curl -fsSL https://bootstrap.pypa.io/get-pip.py -o /opt/get-pip.py \
&& python /opt/get-pip.py \
&& rm -rf /var/lib/apt/lists/* /opt/get-pip.py

# Setup ENV variables
ENV MULTIQC_VERSION="v1.0"

# Install MultiQC
RUN pip install git+git://github.com/ewels/MultiQC.git@$MULTIQC_VERSION

# Create UPPMAX directories
RUN mkdir /pica /proj /sw
