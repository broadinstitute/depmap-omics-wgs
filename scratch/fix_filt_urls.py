import tomllib

import pandas as pd
from google.cloud import storage
from nebelung.terra_workspace import TerraWorkspace

from depmap_omics_wgs.gcs import rewrite_blob

df = pd.read_csv("~/Desktop/filt_urls.tsv", sep="\t")

df["mut_filtered_vcf_fixed"] = (
    df["mut_filtered_vcf"].str.replace(r"filtered.vcf.gz", "")
    + df["entity:sample_id"]
    + ".filtered.vcf.gz"
)

df["mut_filtered_vcf_idx_fixed"] = (
    df["mut_filtered_vcf_idx"].str.replace(r"filtered.vcf.gz.tbi", "")
    + df["entity:sample_id"]
    + ".filtered.vcf.gz.tbi"
)

df["mut_filtering_stats_fixed"] = (
    df["mut_filtering_stats"].str.replace(r"filtering.stats", "")
    + df["entity:sample_id"]
    + ".filtering.stats"
)

storage_client = storage.Client(project="depmap-omics")

bucket = storage.Bucket.from_uri(df["mut_filtered_vcf"].values[0])

for _, r in df.iterrows():
    # rewrite_blob(
    #     src_blob=storage.Blob.from_string(r["mut_filtered_vcf"], client=storage_client),
    #     dest_blob=storage.Blob.from_string(
    #         r["mut_filtered_vcf_fixed"], client=storage_client
    #     ),
    # )

    rewrite_blob(
        src_blob=storage.Blob.from_string(
            r["mut_filtered_vcf_idx"], client=storage_client
        ),
        dest_blob=storage.Blob.from_string(
            r["mut_filtered_vcf_idx_fixed"], client=storage_client
        ),
    )

    rewrite_blob(
        src_blob=storage.Blob.from_string(
            r["mut_filtering_stats"], client=storage_client
        ),
        dest_blob=storage.Blob.from_string(
            r["mut_filtering_stats_fixed"], client=storage_client
        ),
    )

config = {}

with open("./config.toml", "rb") as f:
    config.update(tomllib.load(f))

workspace = TerraWorkspace(
    workspace_namespace=config["terra"]["workspace_namespace"],
    workspace_name=config["terra"]["workspace_name"],
    owners=[],
)

df = df.drop(
    columns=["mut_filtered_vcf", "mut_filtered_vcf_idx", "mut_filtering_stats"]
).rename(
    columns={
        "mut_filtered_vcf_fixed": "mut_filtered_vcf",
        "mut_filtered_vcf_idx_fixed": "mut_filtered_vcf_idx",
        "mut_filtering_stats_fixed": "mut_filtering_stats",
    }
)

workspace.upload_entities(df)
