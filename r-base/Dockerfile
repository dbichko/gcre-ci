FROM quay.io/gcre/base:latest
MAINTAINER "Dmitri Bichko" dbichko@gmail.com

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update -y \
 && apt-get install -y --no-install-suggests --no-install-recommends \
    binutils zlib1g-dev libssl-dev libcurl4-openssl-dev libxml2-dev \
    ca-certificates git r-base-core r-cran-littler \
 && rm -rf /var/lib/apt/lists/*

ENV R_PROFILE /usr/lib/R/etc/Rprofile.site
COPY Rprofile $R_PROFILE

CMD ["bash"]