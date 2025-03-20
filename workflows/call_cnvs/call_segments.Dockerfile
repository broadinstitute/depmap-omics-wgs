FROM bitnami/minideb:bullseye

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get -y update && \
    apt-get -y install --no-install-recommends --no-install-suggests \
    build-essential \
    curl \
    g++ \
    gcc \
    make \
    r-base

RUN R -e "install.packages(c('readr', 'BiocManager'), repos='https://cloud.r-project.org')"
RUN R -e "BiocManager::install('HMMcopy', update=FALSE)"

RUN curl -SL "https://github.com/arq5x/bedtools2/releases/download/v2.31.0/bedtools.static" \
    -o /usr/bin/bedtools && \
    chmod +x /usr/bin/bedtools

COPY ./segment.R /usr/src/
