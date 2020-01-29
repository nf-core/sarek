FROM nfcore/base:1.8
LABEL authors="Maxime Garcia, Szilveszter Juhos" \
      description="Docker image containing all requirements for nf-core/sarek pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
RUN conda env export --name nf-core-sarek-dev > nf-core-sarek-dev.yml
ENV PATH /opt/conda/envs/nf-core-sarek-dev/bin:$PATH