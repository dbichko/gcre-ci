FROM dbichko/gcre:rcpp-clang

LABEL description "benchmark template for geneticsCRE"

ARG arch=native

COPY flags.${arch} /tmp/flags
RUN cat /tmp/flags >> /home/docker/.R/Makevars

RUN mkdir -p install && mv *.out install

RUN git clone https://github.com/dbichko/gcre-ci.git

CMD ["bash"]
