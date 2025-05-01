# depmap-omics-wgs

This repository contains WDL workflows and a Python module to continuously run workflows on new WGS samples.

# Installation

1. Install the required system dependencies:

    - [pyenv](https://github.com/pyenv/pyenv)
    - [Poetry](https://python-poetry.org/)
    - [pre-commit](https://pre-commit.com/)

2. Install the required Python version (3.12.1):

   ```bash
   pyenv install "$(cat .python-version)"
   ```

3. Confirm that `python` maps to the correct version:

   ```
   python --version
   ```

4. Set the Poetry interpreter and install the Python dependencies:

   ```bash
   poetry env use "$(pyenv which python)"
   poetry install
   ```

5. Copy `env.dist` into a new file `.env` and fill it out:

    - `HASURA_URL` can be copied from `ariadne-prod.yml`
    - `HASURA_ADMIN_SECRET` can be copied from the [hasura-admin-secret-prod](https://console.cloud.google.com/security/secret-manager/secret/hasura-admin-secret-prod/versions?project=depmap-gumbo) secret

## Credentials

Download a service account key for [omics-pipeline-runner](https://console.cloud.google.com/iam-admin/serviceaccounts/details/104623708380190564321?project=depmap-omics) and set the `GOOGLE_APPLICATION_CREDENTIALS` environment variable to its location on your filesystem. This will simulate the permissions available inside the remote execution context.

# Execution

TODO

# Development

Run `pre-commit run --all-files` to automatically format your code with [Ruff](https://docs.astral.sh/ruff/) and check static types with [Pyright](https://microsoft.github.io/pyright).

Whenever possible, function/method arguments and return values should be validated with Pydantic or Pandera (if a data frame).

## GraphQL code generation

This repo uses [ariadne-codegen](https://github.com/mirumee/ariadne-codegen) to generate the `gumbo_gql_client` module. It uses the folder of GraphQL queries (`./gql`) and the current GraphQL schema for a particular Gumbo environment to generate all of the Python classes, Pydantic models, and query/mutation methods for interacting with the [Gumbo GraphQL Service](https://github.com/broadinstitute/gumbo_client/tree/main/gumbo-gql-service). To regenerate the module using the current production schema:

```shell
HASURA_ADMIN_SECRET=... poetry run ariadne-codegen --config ariadne-prod.toml
```

# Scratch files

Some Python files in `./scratch` are available to seed the Gumbo task results table.
