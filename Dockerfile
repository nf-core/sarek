FROM nfcore/base:1.13
LABEL authors="Maxime Garcia, Szilveszter Juhos" \
      description="Docker image containing all software requirements for the nf-core/sarek pipeline"

# Install the conda environment
COPY environment.yml /
RUN conda env create --quiet -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-sarek-3.0dev/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-sarek-3.0dev > nf-core-sarek-3.0dev.yml
