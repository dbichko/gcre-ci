FROM dbichko/base

LABEL description "R installation, including common dependencies that are available as pre-built binaries"

USER root

RUN apt-get update -y \
 && apt-get install -y --no-install-suggests --no-install-recommends \
    ca-certificates binutils make git \
    libc-dev zlib1g-dev libssl-dev libcurl4-openssl-dev libxml2-dev \
    r-base-core r-cran-littler \
    r-cran-curl r-cran-openssl \
    r-cran-rcpp r-cran-matrix r-cran-plyr r-cran-dplyr \
    r-cran-stringr r-cran-knitr r-cran-httr r-cran-memoise r-cran-digest r-cran-jsonlite r-cran-withr \
 && rm -rf /var/lib/apt/lists/*

USER docker

ENV R_PROFILE /usr/lib/R/etc/Rprofile.site
COPY Rprofile $R_PROFILE

CMD ["bash"]