# syntax=docker/dockerfile:1

FROM bioconductor/bioconductor_docker:RELEASE_3_16
RUN mkdir /root/rpkg
RUN mkdir /root/examples
COPY . /root/rpkg

RUN OPENAI_API_KEY=${{OPENAI_API_KEY}} R -e "stopifnot(Sys.getenv('OPENAI_API_KEY') != "")"
RUN R -e "devtools::install('/root/rpkg', dependencies=T)"
RUN R -e "devtools::install_github('davisLaboratory/vissE')"
RUN cd /root/rpkg && R -e "library(vissEServer);lapply(list.files('data-raw', '*_example.R', full.names = T), source)"