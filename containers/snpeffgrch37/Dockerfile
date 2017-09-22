FROM maxulysse/snpeff:1.1

LABEL author="Maxime Garcia" \
description="snpEff 4.3m image with GRCH37.75 genome for use in CAW" \
maintainer="maxime.garcia@scilifelab.se"

# Setup ENV variables
ENV GENOME="GRCh37.75"

# Download Genome
RUN java -jar $SNPEFF_HOME/snpEff.jar download $GENOME
