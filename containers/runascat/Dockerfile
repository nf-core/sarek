FROM r-base:3.3.2

LABEL author="Maxime Garcia" \
description="R 3.3.2 image for use in CAW" \
maintainer="maxime.garcia@scilifelab.se"

# Install libraries
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile \
&& R --no-save "install.packages('RColorBrewer')"

# Create UPPMAX directories
RUN mkdir /pica /proj /sw
