FROM dbichko/r-base:latest
MAINTAINER "Dmitri Bichko" dbichko@gmail.com

# some R packages ignore CC and hard-code gcc

RUN apt-get update -y \
 && apt-get install -y --no-install-suggests --no-install-recommends \
    r-cran-curl r-cran-openssl r-cran-httr r-cran-memoise r-cran-digest r-cran-jsonlite r-cran-withr \
    r-cran-rcpp r-cran-matrix r-cran-plyr r-cran-dplyr \
    clang \
 && rm -rf /var/lib/apt/lists/* \
 && ln -s /usr/bin/clang /usr/bin/gcc

COPY Makevars /home/docker/.R/Makevars
RUN chown -R docker:docker /home/docker

USER docker
WORKDIR /home/docker

RUN Rscript -e 'install.packages(c("devtools"), type="source", clean=TRUE, Ncpus=4, verbose=TRUE, quiet=TRUE, keep_outputs=TRUE)' \
 && rm -rf /tmp/downloaded_packages/ /tmp/*/downloaded_packages/ /tmp/*.rds

 # && Rscript -e 'devtools::install_github("hadley/devtools", quiet=TRUE)' \

CMD ["bash"]
