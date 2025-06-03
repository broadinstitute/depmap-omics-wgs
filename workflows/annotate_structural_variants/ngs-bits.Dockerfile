# rewritten from https://github.com/MHH-Humangenetik/ngs-bits-docker/blob/main/Dockerfile
FROM bitnami/minideb:bullseye

# Set NGS-BITS version, default is 'master'. Override with --build-arg NGS_BITS_VERSION=label
ARG NGS_BITS_VERSION=master

# Set environment variables
ENV PATH=/opt/ngs-bits/bin:/bin:$PATH \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8

# Install runtime and build dependencies
RUN install_packages \
    git \
    make \
    g++ \
    qtbase5-dev \
    libqt5xmlpatterns5-dev \
    libqt5sql5-mysql \
    libqt5sql5-odbc \
    libqt5charts5-dev \
    libqt5svg5-dev \
    python3 \
    python3-matplotlib \
    libbz2-dev \
    liblzma-dev \
    libxml2-dev \
    libcurl4 \
    libcurl4-openssl-dev \
    zlib1g-dev \
    pkg-config \
    ca-certificates \
    wget \
    unzip

# Clone and build NGS-Bits
WORKDIR /opt
RUN git clone https://github.com/imgag/ngs-bits.git
WORKDIR /opt/ngs-bits
RUN git checkout ${NGS_BITS_VERSION} && \
    git submodule update --recursive --init && \
    make build_3rdparty && \
    make build_libs_release && \
    make build_tools_release

# Remove build artifacts and only keep built binaries
RUN find /opt/ngs-bits -mindepth 1 -maxdepth 1 ! -name 'bin' -exec rm -rf {} +

# Remove unnecessary packages and install runtime libraries
RUN apt-get remove -y --purge \
        git \
        make \
        g++ \
        qtbase5-dev \
        libqt5xmlpatterns5-dev \
        libqt5charts5-dev \
        libqt5svg5-dev \
        libbz2-dev \
        liblzma-dev \
        libxml2-dev \
        libcurl4-openssl-dev \
        zlib1g-dev && \
    apt-get autoremove -y && \
    install_packages \
        libqt5network5 \
        libqt5xml5 && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
