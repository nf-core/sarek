FROM nfcore/base:1.14
LABEL \
    author="Maxime Garcia" \
    description="VEP image for nf-core pipelines" \
    maintainer="maxime.garcia@scilifelab.se"

# Install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

# Setup default ARG variables
ARG GENOME=GRCh38
ARG SPECIES=homo_sapiens
ARG VEP_VERSION=105
ARG VEP_TAG=105.0

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-vep-${VEP_TAG}/bin:$PATH

# Download Genome
RUN vep_install \
    -a c \
    -c .vep \
    -s ${SPECIES} \
    -y ${GENOME} \
    --CACHE_VERSION ${VEP_VERSION} \
    --CONVERT \
    --NO_BIOPERL --NO_HTSLIB --NO_TEST --NO_UPDATE

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-vep-${VEP_TAG} > nf-core-vep-${VEP_TAG}.yml
