import tomllib

import pandas as pd
from firecloud_api_cds import api as firecloud_api
from google.cloud import storage
from nebelung.terra_workspace import TerraWorkspace
from nebelung.utils import batch_evenly, call_firecloud_api
from tqdm import tqdm

pd.set_option("display.max_columns", 30)
pd.set_option("display.max_colwidth", 50)
pd.set_option("display.max_info_columns", 30)
pd.set_option("display.max_info_rows", 20)
pd.set_option("display.max_rows", 20)
pd.set_option("display.max_seq_items", None)
pd.set_option("display.width", 200)
pd.set_option("expand_frame_repr", True)
pd.set_option("mode.chained_assignment", "warn")


def list_blobs(
    bucket_name: str, prefix: str | None = None, glob: str | None = None
) -> pd.DataFrame:
    """
    Get the names and sizes of existing blobs in a GCS bucket.

    :param bucket_name: the name of the GCS bucket
    :param prefix: an optional prefix for listing
    :param glob: an optional glob for listing
    :return: a data frame of object names and sizes
    """

    storage_client = storage.Client()

    if prefix is not None and glob is not None:
        raise ValueError("At most one of `prefix` and `glob` can be specified")
    elif prefix is not None:
        pages = storage_client.list_blobs(
            bucket_or_name=bucket_name,
            prefix=prefix,
            delimiter="/",
            fields="items(name,crc32c,size,updated),nextPageToken",
        ).pages
    elif glob is not None:
        pages = storage_client.list_blobs(
            bucket_or_name=bucket_name,
            match_glob=glob,
            fields="items(name,crc32c,size,updated),nextPageToken",
        ).pages
    else:
        pages = storage_client.list_blobs(
            bucket_or_name=bucket_name,
            fields="items(name,crc32c,size,updated),nextPageToken",
        ).pages

    blobs = []

    for page in pages:
        blobs.extend(
            [
                {
                    "url": "gs://" + bucket_name + "/" + x.name,
                    "crc32c": x.crc32c,
                    "size": x.size,
                    "gcs_obj_updated_at": x.updated,
                }
                for x in page
            ]
        )

    return pd.DataFrame(blobs)


config = {}

with open("./config.toml", "rb") as f:
    config.update(tomllib.load(f))

legacy_workspace = TerraWorkspace(
    workspace_namespace=config["terra"]["legacy_workspace_namespace"],
    workspace_name=config["terra"]["legacy_workspace_name"],
    owners=[],
)

workspace = TerraWorkspace(
    workspace_namespace=config["terra"]["workspace_namespace"],
    workspace_name=config["terra"]["workspace_name"],
    owners=[],
)

src_samples = legacy_workspace.get_entities("sample")[
    [
        "sample_id",
        "msisensor2_output",
        "msisensor2_output_dis",
        "msisensor2_output_somatic",
        "msisensor2_score",
    ]
]

dest_samples = workspace.get_entities("sample")

src_samples = src_samples.loc[src_samples["sample_id"].isin(dest_samples["sample_id"])]

src_bucket_name = call_firecloud_api(
    firecloud_api.get_workspace,
    namespace=legacy_workspace.workspace_namespace,
    workspace=legacy_workspace.workspace_name,
    fields=["workspace.bucketName"],
)["workspace"]["bucketName"]

dest_bucket_name = call_firecloud_api(
    firecloud_api.get_workspace,
    namespace=workspace.workspace_namespace,
    workspace=workspace.workspace_name,
    fields=["workspace.bucketName"],
)["workspace"]["bucketName"]

storage_client = storage.Client()
src_bucket = storage_client.bucket(src_bucket_name)
dest_bucket = storage_client.bucket(dest_bucket_name)

melted = src_samples.melt(
    id_vars="sample_id", var_name="col", value_name="url"
).dropna()

melted = melted.loc[
    melted["url"].astype("string").str.startswith(f"gs://{src_bucket_name}/")
]

melted["new_url"] = melted["url"].str.replace(
    src_bucket_name, dest_bucket_name, regex=False
)

existing_blobs = list_blobs(dest_bucket_name)
melted_todo = melted.loc[~melted["new_url"].isin(existing_blobs["url"])]

for batch in tqdm(batch_evenly(melted_todo.sort_index(ascending=False), 100)):
    with storage_client.batch(raise_exception=False):
        for _, r in batch.iterrows():
            source_blob = storage.Blob.from_string(r["url"], client=storage_client)
            destination_blob = storage.Blob.from_string(
                r["new_url"], client=storage_client
            )

            blob_copy = src_bucket.copy_blob(
                source_blob, dest_bucket, destination_blob.name
            )

unmelted = (
    melted.pivot(index="sample_id", columns="col", values="new_url")
    .reset_index()
    .astype("string")
)

unmelted = unmelted.merge(
    src_samples[["sample_id", "msisensor2_score"]], how="left", on="sample_id"
)

workspace.upload_entities(unmelted)
