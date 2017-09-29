FROM maxulysse/samtools:1.1

LABEL author="Maxime Garcia" \
description="alleleCount 2.2.0 and SAMTools 1.4 image for use in CAW" \
maintainer="maxime.garcia@scilifelab.se"

# Setup ENV variables
ENV ALLELECOUNT_BIN="v2.2.0.tar.gz" \
ALLELECOUNT_VERSION="2.2.0" \
PATH="$PATH:/opt/bin" \
PERL5LIB="$PERL5LIB:/opt/lib/perl5"

# Install alleleCount
RUN curl -fsSL https://github.com/cancerit/alleleCount/archive/$ALLELECOUNT_BIN -o /opt/$ALLELECOUNT_BIN \
&& tar xvzf /opt/$ALLELECOUNT_BIN -C /opt/ \
&& cd /opt/alleleCount-$ALLELECOUNT_VERSION \
&& ./setup.sh /opt/ \
&& rm /opt/$ALLELECOUNT_BIN
