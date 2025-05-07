FROM bitnami/minideb:bullseye AS python_builder

ENV DEBIAN_FRONTEND="noninteractive" \
    PYTHON_VERSION="3.12.3" \
    PYTHON_UNBUFFERED=1

RUN install_packages \
    curl \
    ca-certificates \
    wget \
    tar \
    gcc \
    make \
    libssl-dev \
    zlib1g-dev \
    libbz2-dev \
    libreadline-dev \
    libsqlite3-dev \
    libffi-dev \
    liblzma-dev

# download and install Python 3.12.3 (prebuilt source tarball from python.org)
RUN curl -fsSL https://www.python.org/ftp/python/${PYTHON_VERSION}/Python-${PYTHON_VERSION}.tgz -o python.tgz && \
    tar -xzf python.tgz && \
    cd Python-${PYTHON_VERSION} && \
    ./configure --prefix=/python --enable-optimizations --with-ensurepip=install && \
    make -j"$(nproc)" && make install && \
    cd / && rm -rf Python-${PYTHON_VERSION} python.tgz

# install open-cravat into a virtual environment
RUN /python/bin/python3 -m venv /venv && \
    /venv/bin/pip install --upgrade pip && \
    /venv/bin/pip install open-cravat==2.12.0


FROM bitnami/minideb:bullseye

ENV DEBIAN_FRONTEND="noninteractive" \
    PATH="/venv/bin:$PATH" \
    PYTHONUNBUFFERED=1 \
    BCFTOOLS_VERSION="1.20"

COPY --from=python_builder /python /python
COPY --from=python_builder /venv /venv

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
    && apt-get -y install --no-install-recommends --no-install-suggests \
        autoconf \
        gcc \
        google-cloud-cli \
        lbzip2 \
        libbz2-dev \
        libcurl4-gnutls-dev \
        libgsl-dev \
        liblzma-dev \
        libperl-dev \
        libssl-dev \
        libz-dev \
        make \
        perl \
        pkg-config \
    && apt-get -y autoremove \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# install bcftools
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
