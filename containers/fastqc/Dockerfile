FROM openjdk:8

LABEL author="Maxime Garcia" \
description="FastQC 0.11.5 image for use in CAW" \
maintainer="maxime.garcia@scilifelab.se"

# Setup ENV variables
ENV FASTQC_BIN="fastqc_v0.11.5.zip"

# Install FastQC
RUN curl -fsSL http://www.bioinformatics.babraham.ac.uk/projects/fastqc/$FASTQC_BIN -o /opt/$FASTQC_BIN \
&& unzip /opt/$FASTQC_BIN -d /opt/ \
&& chmod 755 /opt/FastQC/fastqc \
&& ln -s /opt/FastQC/fastqc /usr/local/bin/fastqc \
&& rm /opt/$FASTQC_BIN

# Create UPPMAX directories
RUN mkdir /pica /proj /sw
