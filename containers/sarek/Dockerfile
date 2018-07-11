FROM nfcore/base:latest

LABEL \
  authors="Maxime.Gracia@scilifelab.se, Szilveszter.Juhos@scilifelab.se" \
  description="Image with tools used in Sarek" \
  maintainers="Maxime.Gracia@scilifelab.se, Szilveszter.Juhos@scilifelab.se"

COPY environment.yml /
RUN conda env update -n root -f /environment.yml && conda clean -a
