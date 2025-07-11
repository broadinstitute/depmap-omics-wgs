FROM python:3.12.9-bullseye

ENV PYTHONUNBUFFERED=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=on \
    PIP_DEFAULT_TIMEOUT=100 \
    PYTHONPATH=/app \
    UV_CACHE_DIR=/tmp/uv-cache

WORKDIR /app

RUN pip install uv

COPY pyproject.toml uv.lock README.md ./
RUN uv sync --no-dev

COPY select_structural_variants ./select_structural_variants
RUN uv sync --no-dev
