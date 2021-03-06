FROM debian:unstable-slim

LABEL maintainer "dbichko@gmail.com"
LABEL description "base image for gcre"

ENV DEBIAN_FRONTEND noninteractive

# some R packages ignore CC and hard-code gcc
RUN apt-get update -y && apt-get dist-upgrade -y \
 && apt-get install -y --no-install-suggests --no-install-recommends \
    apt-utils locales sudo \
    libc-dev zlib1g-dev libssl-dev libcurl4-openssl-dev libxml2-dev \
    ca-certificates binutils make git \
    htop less time nano file gzip bzip2 zip unzip pigz pv wget curl \
    clang libomp-dev libtool \
 && ln -s /usr/bin/clang++ /usr/bin/g++ \
 && rm -rf /var/lib/apt/lists/*

RUN useradd -m docker -s /bin/bash \
 && addgroup docker staff \
 && echo "docker ALL=(ALL:ALL) NOPASSWD:ALL" > /etc/sudoers.d/docker \
 && chmod 0440 /etc/sudoers.d/docker

# envs are inherited in child dockerfiles
ENV TZ America/New_York
ENV LC_ALL en_US.UTF-8
ENV LANG $LC_ALL

RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN echo "$LC_ALL UTF-8" > /etc/locale.gen \
 && mkdir -p /usr/share/locale/ \
 && touch /usr/share/locale/locale.alias \
 && locale-gen $LC_ALL \
 && /usr/sbin/update-locale

USER docker
WORKDIR /home/docker

CMD ["bash"]
