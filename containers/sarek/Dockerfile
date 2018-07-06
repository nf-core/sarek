FROM nfcore/base:latest

LABEL \
  author="Szilveszter Juhos" \
  description="Image with tools used in Sarek" \
  maintainer="szilveszter.juhos@scilifelab.se"

COPY environment.yml /

RUN \
  conda env create -f /environment.yml && \
  conda clean -a

  # Export PATH
ENV PATH /opt/conda/envs/sarek-core/bin:$PATH
