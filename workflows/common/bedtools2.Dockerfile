FROM bitnami/minideb:bullseye

ENV DEBIAN_FRONTEND noninteractive

RUN install_packages bedtools
