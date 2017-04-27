FROM gcre/rcpp

LABEL description "build and install geneticsCRE"

ARG branch
ARG arch

ENV branch ${branch:-master}
ENV arch ${arch:-native}

VOLUME /work

COPY bashrc .bashrc
COPY install .

CMD bootstrap=true bash
