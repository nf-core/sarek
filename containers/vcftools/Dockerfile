FROM nfcore/base:latest

LABEL \
  author="Maxime Garcia" \
  description="vcftools image used in Sarek 2.0" \
  maintainer="maxime.garcia@scilifelab.se"

COPY environment.yml /

RUN \
  conda env create -f /environment.yml && \
  conda clean -a

  # Export PATH
ENV PATH /opt/conda/envs/sarek-vcftools-2.0/bin:$PATH
