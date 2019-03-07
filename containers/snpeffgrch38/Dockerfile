FROM nfcore/base:latest

LABEL \
  author="Maxime Garcia" \
  description="snpEff image with GRCh38.86 genome for use in Sarek" \
  maintainer="maxime.garcia@scilifelab.se"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/sarek-snpeff-2.3/bin:$PATH

# Setup ENV variables
ENV GENOME="GRCh38.86"

# Download Genome
RUN snpEff download -v $GENOME
