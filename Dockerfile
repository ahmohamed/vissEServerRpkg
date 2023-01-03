# syntax=docker/dockerfile:1

FROM bioconductor/bioconductor_docker:RELEASE_3_16
RUN mkdir /root/rpkg
RUN mkdir /root/examples
COPY . /root/rpkg

RUN R -e "devtools::install('/root/rpkg', dependencies=T)"
RUN cd /root/rpkg && R -e "library(vissEServer);lapply(list.files('data-raw', '*_example.R', full.names = T), source)"