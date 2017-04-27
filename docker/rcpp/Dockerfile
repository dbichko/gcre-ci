FROM gcre/r-base

LABEL description "rcpp environment; source packages are built with the chosen compiler (default: clang)"

# only default for now
ARG compiler=clang

USER root

# some R packages ignore CC and hard-code gcc
RUN apt-get update -y \
 && apt-get install -y --no-install-suggests --no-install-recommends \
    clang libomp-dev libtool \
 && ln -s /usr/bin/clang++ /usr/bin/g++ \
 && rm -rf /var/lib/apt/lists/*

COPY Makevars-${compiler} /home/docker/.R/Makevars
COPY flags.* .R/

RUN Rscript -e 'install.packages(c("devtools", "roxygen2"), type="source", clean=TRUE, Ncpus=4, verbose=TRUE, quiet=TRUE, keep_outputs=TRUE)' \
 && rm -rf /tmp/downloaded_packages/ /tmp/*/downloaded_packages/ /tmp/*.rds \
 && cp .R/Makevars .R/Makevars.base \
 && mkdir -p /var/log/rcpp && mv *.out /var/log/rcpp \
 && chown -R docker:docker /home/docker

USER docker

CMD ["bash"]
