FROM python:3.10.19-slim

ENV PYTHONUNBUFFERED=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=on \
    PIP_DEFAULT_TIMEOUT=100 \
    PYTHONPATH=/SignatureAnalyzer-uv \
    UV_VENV=/SignatureAnalyzer-uv/.venv

ENV PATH="$UV_VENV/bin:$PATH"

# install git
RUN apt-get update \
    && apt-get -y install --no-install-recommends --no-install-suggests git \
    && apt-get -y autoremove \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN pip install uv

# clone and install SignatureAnalyzer-uv (forked from SignatureAnalyzer commit
# f48c4a969aa4b65b0833094209ce4340ec3a5375 to get a compatible git submodule for
# SignatureAnalyzer-GPU, which is instead copied in as code, and updated to use uv for
# dependency management)
RUN git clone --recursive https://github.com/broadinstitute/SignatureAnalyzer-uv.git \
    && cd SignatureAnalyzer-uv \
    && uv sync --no-cache

COPY ./run_signature_analyzer.py /app/run_signature_analyzer.py
