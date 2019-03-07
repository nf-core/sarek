FROM nfcore/base:latest

LABEL \
  authors="maxime.garcia@scilifelab.se, szilveszter.juhos@scilifelab.se" \
  description="Image with tools used in Sarek" \
	maintainer="Maxime Garcia <maxime.garcia@scilifelab.se>, Szilveszter Juhos <szilveszter.juhos@scilifelab.se>"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/sarek-2.3/bin:$PATH
