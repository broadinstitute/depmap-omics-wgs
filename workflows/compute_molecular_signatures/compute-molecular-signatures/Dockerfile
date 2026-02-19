FROM python:3.14.2

ENV PYTHONUNBUFFERED=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=on \
    PIP_DEFAULT_TIMEOUT=100 \
    PYTHONPATH=/app \
    UV_VENV=/app/.venv

ENV PATH="$UV_VENV/bin:$PATH"

WORKDIR /app

RUN pip install uv

COPY pyproject.toml uv.lock README.md ./
COPY compute_molecular_signatures/__init__.py ./compute_molecular_signatures/__init__.py
RUN uv sync --no-cache --no-group=dev

COPY compute_molecular_signatures ./compute_molecular_signatures
