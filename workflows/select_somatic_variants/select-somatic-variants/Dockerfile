FROM python:3.12.3-bullseye

ENV PYTHONUNBUFFERED=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=on \
    PIP_DEFAULT_TIMEOUT=100 \
    POETRY_VERSION=1.8.5 \
    POETRY_HOME="/opt/poetry" \
    POETRY_NO_INTERACTION=1 \
    POETRY_VIRTUALENVS_CREATE=true \
    POETRY_VIRTUALENVS_IN_PROJECT=true \
    PYTHONPATH=/app \
    VIRTUAL_ENVIRONMENT_PATH="/app/.venv"

ENV PATH="$VIRTUAL_ENVIRONMENT_PATH/bin:$PATH"

WORKDIR /${PYTHONPATH}

RUN pip install poetry==${POETRY_VERSION}
COPY pyproject.toml poetry.lock ./
RUN poetry install --only main --no-root --no-cache

COPY select_somatic_variants ./select_somatic_variants
RUN poetry install --only-root
