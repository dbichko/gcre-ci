FROM buildpack-deps:stretch

MAINTAINER "Dmitri Bichko" dbichko@gmail.com

RUN apt-get update && apt-get install -y --no-install-recommends \
    apt-utils binutils pv pigz unzip time wget \
    cmake clang \
    less nano locales \
  && rm -rf /var/lib/apt/lists/*

RUN useradd docker \
  && mkdir /home/docker \
  && chown docker:docker /home/docker \
  && addgroup docker staff

RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
  && /usr/share/locale/locale.alias \
  && locale-gen en_US.UTF8 \
  && /usr/sbin/update-locale LANG=en_US.UTF-8

ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

USER docker
WORKDIR /home/docker

CMD ["bash"]