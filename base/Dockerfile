FROM debian:unstable-slim
MAINTAINER "Dmitri Bichko" dbichko@gmail.com

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update -y && apt-get dist-upgrade -y \
 && apt-get install -y --no-install-suggests --no-install-recommends \
    apt-utils locales nano htop less time file gzip bzip2 zip unzip pigz pv wget curl \
    libc-dev libgcc-6-dev make \
 && rm -rf /var/lib/apt/lists/*

ENV TZ America/New_York
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN useradd -m docker -s /bin/bash && addgroup docker staff

RUN echo 'en_US.UTF-8 UTF-8' > /etc/locale.gen \
 && mkdir -p /usr/share/locale/ \
 && touch /usr/share/locale/locale.alias \
 && locale-gen en_US.UTF-8 \
 && /usr/sbin/update-locale LANG=en_US.UTF-8

# envs are inherited in child dockerfiles
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

CMD ["bash"]