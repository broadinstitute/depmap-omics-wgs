FROM python:3.12.3-bullseye

ENV DEBIAN_FRONTEND noninteractive
ENV BCFTOOLS_VERSION="1.20"

RUN apt-get -y update \
    && apt-get -y dist-upgrade \
    && apt-get -y install --no-install-recommends --no-install-suggests \
        autoconf ca-certificates curl gcc libbz2-dev libcurl4-gnutls-dev \
        libgsl-dev liblzma-dev libperl-dev libssl-dev libz-dev make perl \
        pkg-config lbzip2 \
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

ENV BCFTOOLS_PLUGINS /usr/local/src/bcftools/plugins

ENV PYTHONUNBUFFERED=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=on \
    PIP_DEFAULT_TIMEOUT=100 \
    POETRY_VERSION=1.8.3 \
    POETRY_HOME="/opt/poetry" \
    POETRY_NO_INTERACTION=1 \
    POETRY_VIRTUALENVS_CREATE=true \
    POETRY_VIRTUALENVS_IN_PROJECT=true \
    PYTHONPATH=/app \
    VIRTUAL_ENVIRONMENT_PATH="/app/.venv"

ENV PATH="$VIRTUAL_ENVIRONMENT_PATH/bin:$PATH"

WORKDIR /${PYTHONPATH}

RUN pip install poetry
COPY pyproject.toml poetry.lock ./
RUN poetry install --only main --no-root --no-cache

COPY annotate_mutations_postprocess ./annotate_mutations_postprocess
RUN poetry install --only-root

COPY config.toml ./
