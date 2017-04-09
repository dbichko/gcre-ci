FROM dbichko/base:latest

RUN sudo apt-get update -y \
 && sudo apt-get install -y --no-install-suggests --no-install-recommends \
    ca-certificates binutils git \
    zlib1g-dev libssl-dev libcurl4-openssl-dev libxml2-dev \
    r-base-core r-cran-littler \
    r-cran-curl r-cran-openssl \
    r-cran-httr r-cran-memoise r-cran-digest r-cran-jsonlite r-cran-withr \
 && sudo rm -rf /var/lib/apt/lists/*

#    libc-dev libgcc-6-dev make \



ENV R_PROFILE /usr/lib/R/etc/Rprofile.site
COPY Rprofile $R_PROFILE

CMD ["bash"]