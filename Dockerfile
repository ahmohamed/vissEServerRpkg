# syntax=docker/dockerfile:1

FROM bioconductor/bioconductor_docker:RELEASE_3_16
RUN mkdir /root/rpkg
RUN mkdir /root/examples
COPY . /root/rpkg

RUN --mount=type=secret,id=OPENAI_API_KEY \
  R -e "Sys.setenv(OPENAI_API_KEY=readLines('/run/secrets/OPENAI_API_KEY')[[1]]) if(Sys.getenv('OPENAI_API_KEY') != '') stop('none')"

RUN R -e "devtools::install('/root/rpkg', dependencies=T)"
RUN R -e "devtools::install_github('davisLaboratory/vissE')"
RUN cd /root/rpkg && R -e "library(vissEServer);lapply(list.files('data-raw', '*_example.R', full.names = T), source)"