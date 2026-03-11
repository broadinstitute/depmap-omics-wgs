FROM python:3.12.3-bullseye

ENV PYTHONUNBUFFERED=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=on \
    PIP_DEFAULT_TIMEOUT=100 \
    PYTHONPATH=/app \
    UV_VENV=/app/.venv

ENV PATH="$UV_VENV/bin:$PATH"

WORKDIR /${PYTHONPATH}

RUN pip install uv

COPY pyproject.toml uv.lock README.md ./
COPY select_somatic_variants/__init__.py ./select_somatic_variants/__init__.py
RUN uv sync --no-cache --no-group=dev

COPY select_somatic_variants ./select_somatic_variants
