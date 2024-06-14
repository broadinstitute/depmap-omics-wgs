FROM amazoncorretto:22.0.1-alpine

WORKDIR /tmp

RUN apk add --no-cache bash wget unzip && \
    wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip -O snpEff_latest_core.zip && \
    unzip snpEff_latest_core.zip && \
    mkdir -p /app && \
    cp snpEff/snpEff.jar /app/ && \
    cp snpEff/SnpSift.jar /app/ && \
    cp snpEff/snpEff.config /app/ && \
    rm -rf /tmp/* && \
    apk del wget unzip

WORKDIR /app
