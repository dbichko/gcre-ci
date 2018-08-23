FROM gcre/base

LABEL description "build and install geneticsCRE"

ARG compiler=clang
ARG branch
ARG arch

USER root

RUN apt-get update -y \
 && apt-get install -y --no-install-suggests --no-install-recommends \
    r-base-core r-cran-littler r-cran-curl r-cran-openssl \
    r-cran-rcpp r-cran-matrix r-cran-plyr r-cran-dplyr r-cran-igraph r-cran-ggplot2 \
    r-cran-stringr r-cran-knitr r-cran-httr r-cran-memoise r-cran-digest r-cran-jsonlite r-cran-withr \
 && rm -rf /var/lib/apt/lists/*

COPY Makevars-${compiler} /home/docker/.R/Makevars
COPY flags.* .R/

RUN Rscript -e 'install.packages(c("devtools", "roxygen2"), type="source", clean=TRUE, Ncpus=4, verbose=TRUE, quiet=TRUE, keep_outputs=TRUE)' \
 && rm -rf /tmp/downloaded_packages/ /tmp/*/downloaded_packages/ /tmp/*.rds \
 && cp .R/Makevars .R/Makevars.base \
 && mkdir -p /var/log/rcpp && mv *.out /var/log/rcpp \
 && chown -R docker:docker /home/docker

USER docker

ENV R_PROFILE /usr/lib/R/etc/Rprofile.site
COPY Rprofile $R_PROFILE

ENV branch ${branch:-master}
ENV arch ${arch:-native}

VOLUME /work

COPY bashrc .bashrc
COPY install .

CMD bootstrap=true bash
