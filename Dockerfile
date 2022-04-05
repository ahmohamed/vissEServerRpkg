# syntax=docker/dockerfile:1

FROM bioconductor/bioconductor_docker:RELEASE_3_14
RUN R -e "BiocManager::install(c('SpatialExperiment', 'clusterProfiler', 'vissE', ask=F)"

RUN mkdir /root/rpkg
COPY . /root/rpkg

RUN R -e "devtools::install('/root/rpkg', dependencies=T)"