FROM maxulysse/snpeff:1.1

LABEL author="Maxime Garcia" \
description="snpEff 4.3m image with GRCh38.86 genome for use in CAW" \
maintainer="maxime.garcia@scilifelab.se"

# Setup ENV variables
ENV GENOME="GRCh38.86"

# Download Genome
RUN java -jar $SNPEFF_HOME/snpEff.jar download $GENOME
