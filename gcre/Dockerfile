FROM dbichko/gcre:rcpp-clang

LABEL description "benchmark template for geneticsCRE"

VOLUME /results

COPY flags.* .R/
RUN mkdir -p logs && mv *.out logs/
RUN cp .R/Makevars .R/Makevars.base

RUN git clone https://github.com/dbichko/gcre-ci.git

ENV GCRE_ARCH native
ENV GCRE_BRANCH master-bench

# ENTRYPOINT ["gcre-ci/run_benchmark"]
CMD ["gcre-ci/run_benchmark", "full"]
