FROM nfcore/base:latest

LABEL \
  author="Szilveszter Juhos" \
  description="vcfanno Image for Sarek" \
  maintainer="szilveszter.juhos@scilifelab.se"

COPY environment.yml /

RUN \
  conda env create -f /environment.yml && \
  conda clean -a

  # Export PATH
ENV PATH /opt/conda/envs/sarek-vcfanno-2.0/bin:$PATH
