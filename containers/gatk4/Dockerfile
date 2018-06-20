FROM nfcore/base:latest

LABEL \
  author="Maxime Garcia" \
  description="GATK4 Image for Sarek" \
  maintainer="maxime.garcia@scilifelab.se"

COPY environment.yml /

RUN \
  conda env create -f /environment.yml && \
  conda clean -a

  # Export PATH
ENV PATH /opt/conda/envs/sarek-gatk4-2.0/bin:$PATH
