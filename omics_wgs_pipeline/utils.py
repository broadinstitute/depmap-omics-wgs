import json
import uuid
from typing import Any, Callable, Iterable, Type

import pandas as pd
from click import echo
from google.cloud import storage
from nebelung.utils import batch_evenly, type_data_frame

from omics_wgs_pipeline.types import (
    GcsObject,
    PanderaBaseSchema,
    PydanticBaseModel,
    TypedDataFrame,
)


def df_to_model(
    df: pd.DataFrame, pydantic_schema: Type[PydanticBaseModel]
) -> list[PydanticBaseModel]:
    """
    Convert a Pandas data frame to a Pydantic model.

    :param df: a data frame
    :param pydantic_schema: the Pydantic schema to cast the data frame to
    :return: a Pydantic model
    """

    return [pydantic_schema(**x) for x in df.to_dict(orient="records")]


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

    echo(f"Getting metadata about {len(list(urls))} GCS objects")
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
