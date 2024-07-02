FROM debian:bullseye as htslib_builder

ENV VERSION_HTSLIB="1.20"

RUN apt-get update && \
    apt-get -y install build-essential wget zlib1g-dev liblzma-dev libbz2-dev libcurl3-dev

RUN wget -q https://github.com/samtools/htslib/releases/download/${VERSION_HTSLIB}/htslib-${VERSION_HTSLIB}.tar.bz2 && \
    tar -xjf htslib-${VERSION_HTSLIB}.tar.bz2 && \
    cd htslib-${VERSION_HTSLIB} && \
    make -j4 && \
    make install

FROM python:3.12.3-bullseye

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONUNBUFFERED=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=on \
    PIP_DEFAULT_TIMEOUT=100

COPY --from=htslib_builder /usr/local/bin/bgzip /usr/local/bin/bgzip

RUN pip install open-cravat==2.7.3

# install gcloud CLI
RUN apt-get update && \
    apt-get install -y \
    apt-transport-https \
    ca-certificates \
    curl \
    gnupg && \
    echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg add - && \
    apt-get update && \
    apt-get install -y google-cloud-cli && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
