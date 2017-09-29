FROM openjdk:8-slim

LABEL author="Maxime Garcia" \
description="GATK 3.7 image for use in CAW" \
maintainer="maxime.garcia@scilifelab.se"

# Setup ENV variables
ENV GATK_BIN="GenomeAnalysisTK-3.7.tar.bz2" \
GATK_HOME="/opt/gatk"

# Install GATK
COPY $GATK_BIN $GATK_HOME/$GATK_BIN
RUN cd $GATK_HOME \
&& tar -jxf $GATK_BIN \
&& rm $GATK_BIN

# Create UPPMAX directories
RUN mkdir /pica /proj /sw
