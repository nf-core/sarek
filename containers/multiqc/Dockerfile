FROM python:2.7-slim

LABEL \
  author="Maxime Garcia" \
  description="MultiQC 1.1 image used in CAW" \
  maintainer="maxime.garcia@scilifelab.se"

# Build container
ADD \
  build.sh /usr/bin/
RUN \
  build.sh