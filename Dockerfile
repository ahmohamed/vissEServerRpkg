# syntax=docker/dockerfile:1

FROM bioconductor/bioconductor_docker:RELEASE_3_16
RUN mkdir /root/rpkg
RUN mkdir /root/examples
COPY . /root/rpkg

RUN R -e "devtools::install('/root/rpkg', dependencies=T)"
RUN R -e "devtools::install_github('davisLaboratory/vissE')"
RUN --mount=type=secret,id=OPENAI_API_KEY \
  cd /root/rpkg && R -e "Sys.setenv(OPENAI_API_KEY=readLines('/run/secrets/OPENAI_API_KEY')[[1]]);library(vissEServer);lapply(list.files('data-raw', '*_example.R', full.names = T), source)"