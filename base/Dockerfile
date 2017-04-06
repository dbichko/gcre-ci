
FROM debian:unstable-slim

MAINTAINER "Dmitri Bichko" dbichko@gmail.com

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update -y \
  && apt-get install -y --no-install-suggests --no-install-recommends apt-utils \
  && apt-get dist-upgrade -y \
  && apt-get install -y --no-install-suggests --no-install-recommends \
     nano less time file gzip bzip2 zip unzip pigz pv wget curl \
     libc-dev libgcc-6-dev make \
  && rm -rf /var/lib/apt/lists/*


# RUN apt-get update && apt-get install -y --no-install-recommends \
#     libssl-dev libcurl4-openssl-dev \
#     less git nano vim-tiny locales ca-certificates \
#     r-base r-base-dev \
#   && rm -rf /var/lib/apt/lists/*

RUN useradd -m docker -s /bin/bash && addgroup docker staff

CMD ["bash"]