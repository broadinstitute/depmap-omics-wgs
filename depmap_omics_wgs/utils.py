import json
import logging
import uuid
from pathlib import Path
from typing import Any, Callable, Iterable, Type

import pandas as pd
from google.cloud import secretmanager_v1, storage
from nebelung.terra_workflow import TerraWorkflow
from nebelung.terra_workspace import TerraWorkspace
from nebelung.utils import batch_evenly, type_data_frame

from depmap_omics_wgs.types import (
    GcsObject,
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


def get_gcs_object_metadata(
    urls: Iterable[str], gcp_project_id: str
) -> TypedDataFrame[GcsObject]:
    """
    Check existence and get metadata (size, hash, etc.) of GCS objects.

    :param urls: iterable of GCS URLs
    :param gcp_project_id: the ID of a GCP project to use for billing
    :return: data frame of object URLs and metadata
    """

    logging.info(f"Getting metadata about {len(list(urls))} GCS objects")
    storage_client = storage.Client(project=gcp_project_id)
    blobs = {}

    # the GCS batch context below has a max batch size of 1000, so do this outer layer
    # of batching, too)
    for batch in batch_evenly(urls, max_batch_size=200):
        with storage_client.batch(raise_exception=False):
            for url in batch:
                blob = storage.Blob.from_string(url, client=storage_client)
                bucket = storage_client.bucket(
                    blob.bucket.name, user_project=gcp_project_id
                )
                blob = bucket.get_blob(blob.name)
                blobs[url] = blob

    df = pd.DataFrame(
        [
            {
                "url": k,
                "size": v.size,
                "crc32c_hash": v.crc32c,
                "created_at": v.time_created,
            }
            for k, v in blobs.items()
        ]
    )

    # batching without raising exceptions makes all columns NA if an object was missing
    df = df.dropna()

    return type_data_frame(df, GcsObject)


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
        **kwargs,
    )


def submit_delta_job(
    terra_workspace: TerraWorkspace,
    terra_workflow: TerraWorkflow,
    entity_type: str,
    entity_set_type: str,
    entity_id_col: str,
    expression: str,
    resubmit_n_times: int = 1,
    dry_run: bool = True,
    input_cols: set[str] | None = None,
    output_cols: set[str] | None = None,
):
    """
    Identify entities in a Terra data table that need to have a workflow run on them by:

        1. checking for the presence of a workflow output in a data table column
        2. confirming the entity is eligible to be submitted in a job by checking for
           previous submissions of that same entity to the workflow

    :param terra_workspace: a `TerraWorkspace` instance
    :param terra_workflow: a TerraWorkflow instance for the method
    :param entity_type: the name of the Terra entity type
    :param entity_set_type: the name of the Terra entity set type for `entity_type`
    :param entity_id_col: the name of the ID column for the entity type
    :param expression: the entity type expression (e.g. "this.samples")
    :param resubmit_n_times: the number of times to resubmit an entity in the event it
    has failed in the past
    :param dry_run: whether to skip updates to external data stores
    :param input_cols: the set of column names that must all be present in the entity
    type in order for an entity to be submittable
    :param output_cols: the set of column names that must all be missing in the entity
    type in order for an entity to be submittable
    """

    # get the method config for this workflow in this workspace
    workflow_config = terra_workspace.get_workflow_config(terra_workflow)

    assert not workflow_config["deleted"]
    assert workflow_config["rootEntityType"] == entity_type

    # identify columns in data table used for input/output if not explicitly provided
    if input_cols is None:
        input_cols = {
            v[5:] for k, v in workflow_config["inputs"].items() if v.startswith("this.")
        }

    if output_cols is None:
        output_cols = {
            v[5:]
            for k, v in workflow_config["outputs"].items()
            if v.startswith("this.")
        }

    # get the entities for this workflow entity type
    entities = terra_workspace.get_entities(entity_type)

    # ensure columns exist to check for populated values
    for c in input_cols.union(output_cols):
        if c not in entities.columns:
            entities[c] = pd.NA
        else:
            entities[c] = entities[c].replace({"": pd.NA})

    # identify entities that have all required inputs but no outputs
    entities_todo = entities.loc[
        entities[list(input_cols)].notna().all(axis=1)
        & entities[list(output_cols)].isna().all(axis=1)
    ]

    if len(entities_todo) == 0:
        logging.info(f"No {entity_type}s to run {terra_workflow.method_name} for")
        return

    # get statuses of submitted entity workflow statuses
    submittable_entities = terra_workspace.check_submittable_entities(
        entity_type,
        entity_ids=entities_todo[entity_id_col],
        terra_workflow=terra_workflow,
        resubmit_n_times=resubmit_n_times,
        force_retry=False,
    )

    logging.info(f"Submittable entities: {submittable_entities}")

    if len(submittable_entities["failed"]) > 0:
        raise RuntimeError("Some entities have failed too many times")

    # don't submit jobs for entities that are currently running, completed, or failed
    # too many times
    entities_todo = entities_todo.loc[
        entities_todo[entity_id_col].isin(
            list(
                submittable_entities["unsubmitted"]
                .union(submittable_entities["retryable"])
                .difference(submittable_entities["running"])
            )
        )
    ]

    if dry_run:
        logging.info(f"(skipping) Submitting {terra_workflow.method_name} job")
        return

    entity_set_id = terra_workspace.create_entity_set(
        entity_type,
        entity_ids=entities_todo[entity_id_col],
        suffix=terra_workflow.method_name,
    )

    terra_workspace.submit_workflow_run(
        terra_workflow=terra_workflow,
        entity=entity_set_id,
        etype=entity_set_type,
        expression=expression,
        use_callcache=True,
        use_reference_disks=False,
        memory_retry_multiplier=1.5,
    )
