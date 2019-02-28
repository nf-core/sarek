FROM nfcore/base:latest

LABEL \
  author="Maxime Garcia" \
  description="snpEff image with GRCh37.75 genome for use in Sarek" \
  maintainer="maxime.garcia@scilifelab.se"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/sarek-snpeff-2.3/bin:$PATH

# Setup ENV variables
ENV GENOME="GRCh37.75"

# Download Genome
RUN snpEff download -v $GENOME
