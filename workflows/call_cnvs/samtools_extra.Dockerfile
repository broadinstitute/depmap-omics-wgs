FROM us-central1-docker.pkg.dev/depmap-omics/terra-images/samtools:production

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get -y update

RUN curl -SL "https://github.com/arq5x/bedtools2/releases/download/v2.31.0/bedtools.static" \
    -o /usr/bin/bedtools && \
    chmod +x /usr/bin/bedtools

RUN curl -SL "https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/bedGraphToBigWig" \
    -o /usr/bin/bedGraphToBigWig && \
    chmod +x /usr/bin/bedGraphToBigWig

RUN curl -SL "https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/bigWigToWig" \
    -o /usr/bin/bigWigToWig && \
    chmod +x /usr/bin/bigWigToWig
