FROM bitnami/minideb:bullseye

RUN apt-get update
RUN apt-get -y upgrade
RUN apt-get install -y git libgomp1

WORKDIR /app

RUN git clone https://github.com/niu-lab/msisensor2.git
RUN cd msisensor2
RUN chmod +x msisensor2
ENV PATH="/app/msisensor2:${PATH}"
