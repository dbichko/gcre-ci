
# official r base
FROM rocker/r-base:latest

MAINTAINER "Dmitri Bichko" dbichko@gmail.com

RUN apt-get update && apt-get install -t unstable -y --no-install-recommends \
    apt-utils binutils pv pigz unzip time wget curl \
    build-essential autoconf automake cmake clang gcc \
    libssl-dev libcurl4-openssl-dev libxml2-dev \
    git nano \
  && rm -rf /var/lib/apt/lists/*

RUN Rscript -e "install.packages(c('plyr', 'dplyr', 'devtools'))" \
  && Rscript -e "devtools::install_github('hadley/devtools')" \
  && rm -rf /tmp/downloaded_packages/ /tmp/*.rds

USER docker
WORKDIR /home/docker

CMD ["bash"]