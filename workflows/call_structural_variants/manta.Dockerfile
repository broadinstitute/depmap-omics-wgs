FROM centos:6.10

RUN curl -OL https://github.com/Illumina/manta/releases/download/v1.6.0/manta-1.6.0.centos6_x86_64.tar.bz2
RUN tar xjf manta-1.6.0.centos6_x86_64.tar.bz2
RUN rm manta-1.6.0.centos6_x86_64.tar.bz2

ENV PATH="/manta-1.6.0.centos6_x86_64/bin:${PATH}"
