FROM r-base:3.3.2

LABEL \
  author="Maxime Garcia" \
  description="R 3.3.2 image for use in Sarek" \
  maintainer="maxime.garcia@scilifelab.se"

# Install libraries
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile \
&& Rscript -e "install.packages('RColorBrewer')"
