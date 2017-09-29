FROM ubuntu:16.04

LABEL author="Maxime Garcia" \
description="VEP 90.1 image for use in CAW" \
maintainer="maxime.garcia@scilifelab.se"

# Install libraries
RUN apt-get update && apt-get install -y --no-install-recommends \
apache2 \
build-essential \
ca-certificates \
cpanminus \
curl \
git \
libmysqlclient-dev \
libpng12-dev \
libssl-dev \
manpages \
mysql-client \
openssl \
perl \
perl-base \
unzip \
vim \
wget \
&& rm -rf /var/lib/apt/lists/* \
&& cpanm DBI DBD::mysql

# Set up ENV variables
ENV CFLAGS="-fPIC" \
DEPS="/home/vep/src" \
HOME="/home/vep" \
HTSLIB_DIR="/home/vep/src/htslib" \
KENT_SRC="/home/vep/src/kent-335_base/src" \
MACHTYPE="x86_64" \
PATH="/home/vep/src/ensembl-vep:$PATH" \
PERL5LIB="$PERL5LIB:/home/vep/src/bioperl-live-release-1-6-924:/home/vep/src/ensembl/modules/:/home/vep/src/ensembl-vep/modules/:/home/vep/src/ensembl-variation/modules/:/home/vep/src/ensembl-funcgen/modules/:/home/vep/src/ensembl-io/modules/"

# Preparing install
RUN mkdir -p $DEPS \
&& cd $DEPS \
&& git clone https://github.com/Ensembl/ensembl-funcgen ensembl-funcgen \
&& git clone https://github.com/Ensembl/ensembl-io ensembl-io \
&& git clone https://github.com/Ensembl/ensembl-variation ensembl-variation \
&& git clone https://github.com/Ensembl/ensembl-vep ensembl-vep \
&& git clone https://github.com/Ensembl/ensembl ensembl \
&& $DEPS/ensembl-vep/travisci/get_dependencies.sh \
&& $DEPS/ensembl-vep/travisci/build_c.sh \
&& cd $HTSLIB_DIR \
&& make install \
&& cpanm --installdeps --with-recommends --notest --cpanfile $DEPS/ensembl/cpanfile . \
&& cpanm --installdeps --with-recommends --notest --cpanfile $DEPS/ensembl-vep/cpanfile .

# Installing VEP
RUN chmod u+x $DEPS/ensembl-vep/*.pl
RUN $DEPS/ensembl-vep/INSTALL.pl -a a -l

# Create UPPMAX directories
RUN mkdir /pica /proj /sw
