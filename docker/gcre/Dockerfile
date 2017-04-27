FROM gcre/rcpp

LABEL description "build and install geneticsCRE"

ENV branch master
ENV arch native

VOLUME /work

COPY bashrc .bashrc
COPY install .

CMD bootstrap=true bash
