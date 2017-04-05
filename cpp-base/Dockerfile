
FROM debian:unstable-slim

MAINTAINER "Dmitri Bichko" dbichko@gmail.com

RUN apt-get update && apt-get install -y --no-install-recommends \
    apt-utils binutils pv pigz unzip time wget curl \
    build-essential autoconf automake cmake clang gcc \
    libssl-dev libcurl4-openssl-dev \
    less git nano vim-tiny locales ca-certificates \
    r-base r-base-dev \
  && rm -rf /var/lib/apt/lists/*

RUN useradd docker \
  && mkdir /home/docker \
  && chown docker:docker /home/docker \
  && addgroup docker staff

RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
  && locale-gen en_US.utf8 \
  && /usr/sbin/update-locale LANG=en_US.UTF-8

ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

USER docker
WORKDIR /home/docker

CMD ["bash"]