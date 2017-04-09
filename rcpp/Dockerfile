FROM dbichko/r-base

LABEL description "rcpp environment; source packages are built with the chosen compiler (default: clang)"

ARG compiler=clang

USER root

# some R packages ignore CC and hard-code gcc
RUN apt-get update -y \
 && apt-get install -y --no-install-suggests --no-install-recommends \
    clang libtool \
 && ln -s /usr/bin/clang++ /usr/bin/g++ \
 && rm -rf /var/lib/apt/lists/*

COPY Makevars-${compiler} /home/docker/.R/Makevars
RUN chown -R docker:docker /home/docker

USER docker

RUN Rscript -e 'install.packages(c("devtools", "roxygen2"), type="source", clean=TRUE, Ncpus=4, verbose=TRUE, quiet=TRUE, keep_outputs=TRUE)' \
 && rm -rf /tmp/downloaded_packages/ /tmp/*/downloaded_packages/ /tmp/*.rds

CMD ["bash"]
