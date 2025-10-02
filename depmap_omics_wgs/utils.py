import json
import logging
import uuid
from pathlib import Path
from typing import Any, Callable, Type

import pandas as pd
from google.cloud import secretmanager_v1
from nebelung.terra_workflow import TerraWorkflow
from nebelung.utils import type_data_frame

from depmap_omics_wgs.types import (
    PanderaBaseSchema,
    PydanticBaseModel,
    TypedDataFrame,
)


def get_hasura_creds(gumbo_env: str) -> dict[str, str]:
    """
    Get URL and password for Hasura GraphQL API.

    :param gumbo_env: the Gumbo env to get credentials for ('staging' or 'prod')
    :return: a dictionary with the GraphQL API URL and password
    """

    logging.info(f"Getting Hasura credentials for {gumbo_env}")

    return {
        "url": get_secret_from_sm(
            f"projects/814840278102/secrets/hasura-{gumbo_env}-api-url/versions/latest"
        ),
        "password": get_secret_from_sm(
            f"projects/814840278102/secrets/hasura-admin-secret-{gumbo_env}/versions/latest"
        ),
    }


def get_secret_from_sm(name: str) -> str:
    """
    Get the value of a secret from GCP Secret Manager.

    :param name: the fully-qualified name of a secret
    :return: the secret's decoded value
    """

    client = secretmanager_v1.SecretManagerServiceClient()
    request = secretmanager_v1.AccessSecretVersionRequest(mapping={"name": name})
    response = client.access_secret_version(request=request)
    return response.payload.data.decode()


def model_to_df(
    model: PydanticBaseModel,
    pandera_schema: Type[PanderaBaseSchema],
    records_key: str = "records",
    remove_unknown_cols: bool = False,
    mutator: Callable[[pd.DataFrame], pd.DataFrame] = lambda _: _,
) -> TypedDataFrame[PanderaBaseSchema]:
    """
    Dump a Pydantic model and convert it to a data frame typed by a Pandera schema.

    :param model: a Pydandict model containing a list of objects keyed by `records_key`
    :param pandera_schema: the Pandera schema to cast the model to
    :param records_key: the key/method name in `model` containing the records
    :param remove_unknown_cols: remove columns not specified in the schema
    :param mutator: an optional function to call on the data frame before typing (e.g.
    to rename columns to expected Pydantic field names)
    """

    records = model.model_dump()[records_key]

    df = pd.DataFrame(records)
    df = mutator(df)
    return type_data_frame(df, pandera_schema, remove_unknown_cols)


def compute_uuidv3(
    d: dict[str, Any], uuid_namespace: str, keys: set[str] | None = None
) -> str:
    """
    Compute a consistent UUID-formatted ID for a dictionary.

    :param d: a dictionary
    :param uuid_namespace: a namespace for generated UUIDv3s
    :param keys: an optional subset of keys to use for hashing, otherwise use all keys
    :return: the UUIDv3 as a string
    """

    if keys is None:
        d_subset = d
    else:
        d_subset = {k: v for k, v in d.items() if k in keys}

    hash_key = json.dumps(d_subset, sort_keys=True)
    return str(uuid.uuid3(uuid.UUID(uuid_namespace), hash_key))


def make_workflow_from_config(
    config: dict[str, Any], workflow_name: str, **kwargs: Any
) -> TerraWorkflow:
    """
    Make a TerraWorkflow object from a config entry.

    :param config: a config dictionary
    :param workflow_name: the name of the workflow referenced in the config
    :return: a TerraWorkflow instance
    """

    return TerraWorkflow(
        method_namespace=config["terra"][workflow_name]["method_namespace"],
        method_name=config["terra"][workflow_name]["method_name"],
        method_config_namespace=config["terra"][workflow_name][
            "method_config_namespace"
        ],
        method_config_name=config["terra"][workflow_name]["method_config_name"],
        method_synopsis=config["terra"][workflow_name]["method_synopsis"],
        workflow_wdl_path=Path(
            config["terra"][workflow_name]["workflow_wdl_path"]
        ).resolve(),
        method_config_json_path=Path(
            config["terra"][workflow_name]["method_config_json_path"]
        ).resolve(),
        workflow_inputs_json_path=Path(
            config["terra"][workflow_name]["workflow_inputs_json_path"]
        ).resolve(),
        **kwargs,
    )
