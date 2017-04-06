FROM quay.io/gcre/r-base:latest
MAINTAINER "Dmitri Bichko" dbichko@gmail.com

RUN apt-get update -y \
 && apt-get install -y --no-install-suggests --no-install-recommends gcc libtool \
 && rm -rf /var/lib/apt/lists/*

COPY Makevars /home/docker/.R/Makevars
RUN chown -R docker:docker /home/docker

USER docker
WORKDIR /home/docker

RUN Rscript -e 'install.packages(c("devtools"), type="source", clean=TRUE, Ncpus=4, verbose=TRUE, quiet=TRUE, keep_outputs=TRUE)' \
 && Rscript -e 'devtools::install_github("hadley/devtools", quiet=TRUE)' \
 && rm -rf /tmp/downloaded_packages/ /tmp/*/downloaded_packages/ /tmp/*.rds

CMD ["bash"]
