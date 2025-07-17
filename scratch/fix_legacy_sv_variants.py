import tomllib

import pandas as pd
from firecloud_api_cds import api as firecloud_api
from google.cloud import storage
from nebelung.terra_workspace import TerraWorkspace
from nebelung.utils import call_firecloud_api

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

workspace = TerraWorkspace(
    workspace_namespace=config["terra"]["workspace_namespace"],
    workspace_name=config["terra"]["workspace_name"],
    owners=[],
)

bucket_name = call_firecloud_api(
    firecloud_api.get_workspace,
    namespace=workspace.workspace_namespace,
    workspace=workspace.workspace_name,
    fields=["workspace.bucketName"],
)["workspace"]["bucketName"]

samples = (
    workspace.get_entities("sample")
    .loc[:, ["sample_id", "sv_selected_somatic_sv"]]
    .dropna()
)

samples = samples.loc[samples["sv_selected_somatic_sv"].str.endswith(".bedpe")]

samples["sv_selected_somatic_sv_fixed"] = samples["sv_selected_somatic_sv"].str.replace(
    ".svs.expanded.reannotated.filtered.bedpe", ".selected_somatic_sv.parquet"
)

blobs = list_blobs(bucket_name, glob="**/*.selected_somatic_sv.parquet")
samples["done"] = samples["sv_selected_somatic_sv_fixed"].isin(blobs["url"])
samples = samples.sort_values("sample_id")

min_depth = 5


def fix(r):
    print(r["sample_id"])
    df = pd.read_csv(r["sv_selected_somatic_sv"], sep="\t")

    df["dp"] = df["PR"].str.split(",", expand=True).astype("Int64").sum(axis=1) + df[
        "SR"
    ].str.split(",", expand=True).astype("Int64").sum(axis=1)

    df = df.loc[df["dp"].ge(min_depth)]

    df.to_parquet(r["sv_selected_somatic_sv_fixed"])


for _, r in samples.loc[~samples["done"]].iterrows():
    fix(r)

workspace.upload_entities(
    samples.drop(columns=["sv_selected_somatic_sv", "done"]).rename(
        columns={"sv_selected_somatic_sv_fixed": "sv_selected_somatic_sv"}
    )
)
