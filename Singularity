From:nfcore/base
Bootstrap:docker

%labels
    MAINTAINER Maxime Garcia <maxime.garcia@scilifelab.se>
    DESCRIPTION Singularity image containing all requirements for the Sarek pipeline
    VERSION 2.1.0

%environment
    PATH=/opt/conda/envs/sarek-2.2.1/bin:$PATH
    export PATH

%files
    environment.yml /

%post
    /opt/conda/bin/conda env create -f /environment.yml
    /opt/conda/bin/conda clean -a
