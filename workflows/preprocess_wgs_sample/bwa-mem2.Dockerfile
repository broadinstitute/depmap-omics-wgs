FROM bitnami/minideb:bullseye

ENV DEBIAN_FRONTEND noninteractive
ENV BWA_MEM2_VERSION="2.2.1"
ENV SAMTOOLS_VERSION="1.20"

RUN set -e \
    && apt-get -y update \
    && apt-get -y install --no-install-recommends --no-install-suggests \
    ca-certificates \
    autoconf \
    automake \
    make \
    gcc \
    perl \
    zlib1g-dev \
    libbz2-dev \
    lbzip2 \
    liblzma-dev \
    libcurl4-gnutls-dev \
    libssl-dev \
    libncurses5-dev \
    libdeflate-dev \
    curl

RUN curl -SL \
    https://github.com/bwa-mem2/bwa-mem2/releases/download/v${BWA_MEM2_VERSION}/bwa-mem2-${BWA_MEM2_VERSION}_x64-linux.tar.bz2 \
    -o /tmp/bwa-mem2.tar.bz2 \
    && tar xvf /tmp/bwa-mem2.tar.bz2 -C /usr/local/src --remove-files \
    && mv /usr/local/src/bwa-mem2-* /usr/local/src/bwa-mem2 \
    && find /usr/local/src/bwa-mem2 -maxdepth 1 -type f -executable \
    -exec ln -s {} /usr/local/bin \;

RUN curl -SL \
    https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    -o /tmp/samtools.tar.bz2 \
    && tar xvf /tmp/samtools.tar.bz2 -C /usr/local/src --remove-files \
    && mv /usr/local/src/samtools-* /usr/local/src/samtools \
    && cd /usr/local/src/samtools/htslib-* \
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
