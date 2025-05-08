FROM python:3.12.3-bullseye

ENV DEBIAN_FRONTEND="noninteractive" \
    PYTHONUNBUFFERED=1 \
    PIP_DISABLE_PIP_VERSION_CHECK="on" \
    PIP_DEFAULT_TIMEOUT=100 \
    BCFTOOLS_VERSION="1.20"

RUN pip install open-cravat==2.12.0

# install gcloud CLI
RUN apt-get -y dist-upgrade \
    && apt-get -y update \
    && apt-get -y install --no-install-recommends --no-install-suggests \
        apt-transport-https \
        ca-certificates \
        curl \
        gnupg \
    && curl -fsSL https://packages.cloud.google.com/apt/doc/apt-key.gpg | \
        gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg \
    && echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" \
        > /etc/apt/sources.list.d/google-cloud-sdk.list \
    && apt-get -y update \
    && apt-get -y install --no-install-recommends --no-install-suggests google-cloud-cli

# install bcftools
RUN apt-get -y install --no-install-recommends --no-install-suggests \
        autoconf \
        gcc \
        lbzip2 \
        libbz2-dev \
        libcurl4-gnutls-dev \
        libffi-dev \
        libgsl-dev \
        liblzma-dev \
        libperl-dev \
        libreadline-dev \
        libsqlite3-dev \
        libssl-dev \
        libz-dev \
        make \
        perl \
        pkg-config \
        tar \
        zlib1g-dev \
    && apt-get -y autoremove \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN curl -SL \
    https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 \
    -o /tmp/bcftools.tar.bz2 \
    && tar xvf /tmp/bcftools.tar.bz2 -C /usr/local/src --remove-files \
    && mv /usr/local/src/bcftools-* /usr/local/src/bcftools \
    && cd /usr/local/src/bcftools/htslib-* \
    && autoheader \
    && autoconf \
    && ./configure \
    && make \
    && make install \
    && cd .. \
    && autoheader \
    && autoconf \
    && ./configure --enable-libgsl --enable-perl-filters \
    && make \
    && make install

ENV BCFTOOLS_PLUGINS="/usr/local/src/bcftools/plugins"
