FROM nfcore/base:latest

LABEL \
  authors="Maxime.Garcia@scilifelab.se, Szilveszter.Juhos@scilifelab.se" \
  description="Image with tools used in Sarek" \
	maintainer="Maxime Garcia <maxime.garcia@scilifelab.se>, Szilveszter Juhos <Szilveszter.Juhos@scilifelab.se>"

COPY environment.yml /
RUN conda env update -n root -f /environment.yml && conda clean -a
ENV PATH /opt/conda/bin:$PATH
