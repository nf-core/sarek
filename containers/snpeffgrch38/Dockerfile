FROM maxulysse/snpeff:latest

LABEL \
  author="Maxime Garcia" \
  description="snpEff image with GRCh38.86 genome for use in Sarek" \
  maintainer="maxime.garcia@scilifelab.se"

# Setup ENV variables
ENV GENOME="GRCh38.86"

# Download Genome
RUN java -jar $SNPEFF_HOME/snpEff.jar download $GENOME
