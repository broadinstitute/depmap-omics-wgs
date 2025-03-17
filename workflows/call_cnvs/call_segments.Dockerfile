FROM bitnami/minideb:bullseye

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get -y update && \
    apt-get -y install --no-install-recommends --no-install-suggests \
    r-base

RUN R -e "install.packages(c('readr', 'HMMcopy'), repos='https://cloud.r-project.org', quiet=TRUE)"

COPY ./segment.R /usr/src/
