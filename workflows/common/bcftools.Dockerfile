FROM rust:1.81-bullseye AS gc_prop_builder

RUN apt-get -y update \
    && apt-get -y dist-upgrade \
    && apt-get install -y --no-install-recommends --no-install-suggests libclang-dev

WORKDIR /usr/src

RUN cargo new --bin gc_prop_in_window
WORKDIR /usr/src/gc_prop_in_window

COPY ./gc_prop_in_window/Cargo.toml ./gc_prop_in_window/Cargo.lock ./
RUN cargo build --release
RUN rm -rf src

COPY ./gc_prop_in_window/src ./src
RUN rm ./target/release/deps/gc_prop_in_window*
RUN cargo build --release

FROM bitnami/minideb:bullseye

WORKDIR /app
COPY --from=gc_prop_builder /usr/src/gc_prop_in_window/target/release/gc_prop_in_window .

ENV DEBIAN_FRONTEND="noninteractive"
ENV BCFTOOLS_VERSION="1.20"

ADD https://bootstrap.pypa.io/get-pip.py /tmp/get-pip.py

RUN apt-get -y update \
    && apt-get -y dist-upgrade \
    && apt-get -y install --no-install-recommends --no-install-suggests \
        autoconf ca-certificates curl file gcc libbz2-dev libcurl4-gnutls-dev \
        libgsl-dev liblzma-dev libperl-dev libssl-dev libz-dev make perl \
        pkg-config python3-dev python3-distutils lbzip2 \
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

# also install csvkit
RUN set -e \
    && /usr/bin/python3 /tmp/get-pip.py \
    && pip install -U --no-cache-dir csvkit \
    && rm -f /tmp/get-pip.py

ENV BCFTOOLS_PLUGINS="/usr/local/src/bcftools/plugins"
