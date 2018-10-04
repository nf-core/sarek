FROM nfcore/base:latest

LABEL \
  authors="Maxime.Garcia@scilifelab.se, Szilveszter.Juhos@scilifelab.se" \
  description="Image with tools used in Sarek" \
	maintainer="Maxime Garcia <maxime.garcia@scilifelab.se>, Szilveszter Juhos <Szilveszter.Juhos@scilifelab.se>"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/sarek-2.2.1/bin:$PATH
