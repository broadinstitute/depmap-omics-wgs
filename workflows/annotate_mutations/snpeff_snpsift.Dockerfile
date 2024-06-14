FROM alpine:3.19 as htslib_builder

ENV VERSION_HTSLIB="1.20"

RUN apk update && apk add --no-cache curl curl-dev gcc make libc-dev ncurses-dev zlib-dev xz-dev bzip2-dev

RUN wget -q https://github.com/samtools/htslib/releases/download/${VERSION_HTSLIB}/htslib-${VERSION_HTSLIB}.tar.bz2 && \
    tar -xjf htslib-${VERSION_HTSLIB}.tar.bz2 && \
    cd htslib-${VERSION_HTSLIB} && \
    make -j4 && \
    make install

FROM amazoncorretto:22.0.1-alpine3.19

WORKDIR /tmp

RUN apk update && apk add --no-cache libcurl xz-dev bzip2-dev bash wget unzip && \
    wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip -O snpEff_latest_core.zip && \
    unzip snpEff_latest_core.zip && \
    mkdir -p /app && \
    cp snpEff/snpEff.jar /app/ && \
    cp snpEff/SnpSift.jar /app/ && \
    cp snpEff/snpEff.config /app/ && \
    rm -rf /tmp/* && \
    apk del wget unzip

COPY --from=htslib_builder /usr/local/bin/bgzip /usr/local/bin/bgzip

WORKDIR /app
