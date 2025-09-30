FROM bitnami/minideb:bullseye

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get -y dist-upgrade \
    && apt-get -y update \
    && apt-get -y install --no-install-recommends --no-install-suggests \
        apt-transport-https \
        ca-certificates \
        curl \
        gnupg \
    && echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" \
        | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list \
    && curl https://packages.cloud.google.com/apt/doc/apt-key.gpg \
        | apt-key --keyring /usr/share/keyrings/cloud.google.gpg add - \
    && apt-get -y update \
    && apt-get -y install --no-install-recommends --no-install-suggests \
        autoconf \
        build-essential \
        gcc \
        google-cloud-cli \
        lbzip2 \
        libbz2-dev \
        libcurl4-gnutls-dev \
        libdeflate-dev \
        liblz4-dev \
        liblzma-dev \
        libncurses5-dev \
        libz-dev \
        make \
        xz-utils \
    && apt-get -y autoremove \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

ENV SAMTOOLS_VERSION="1.22.1"

RUN curl -SL https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    -o /tmp/samtools.tar.bz2 \
    && tar xvf /tmp/samtools.tar.bz2 -C /usr/local/src --remove-files \
    && mv /usr/local/src/samtools-* /usr/local/src/samtools \
    && cd /usr/local/src/samtools \
    && autoheader \
    && autoconf -Wno-syntax \
    && ./configure --with-htslib \
    && make \
    && make install
