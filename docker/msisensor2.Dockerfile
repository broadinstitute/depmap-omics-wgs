FROM bitnami/minideb:bullseye

RUN apt-get update
RUN apt-get -y upgrade
RUN apt-get install -y build-essential curl git libbz2-dev libcurl3-dev liblzma-dev \
    libgsl-dev libncurses5-dev wget zip

RUN git clone https://github.com/niu-lab/msisensor2.git
RUN cd msisensor2
RUN chmod +x msisensor2
ENV PATH="/msisensor2:${PATH}"
