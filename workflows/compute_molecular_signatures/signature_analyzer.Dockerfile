FROM python:3.10-slim

# install git
RUN apt-get update \
    && apt-get install -y git \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# upgrade pip and install required Python tools
RUN pip install --upgrade pip setuptools wheel

# clone and install specific commit from SignatureAnalyzer (need to use an old commit
# to get a compatible git submodule for SignatureAnalyzer-GPU)
RUN git clone --recursive https://github.com/broadinstitute/getzlab-SignatureAnalyzer.git \
    && cd getzlab-SignatureAnalyzer \
    && git checkout f48c4a969aa4b65b0833094209ce4340ec3a5375 \
    && pip install -e ./

WORKDIR /app

COPY ./run_signature_analyzer.py /app/run_signature_analyzer.py
