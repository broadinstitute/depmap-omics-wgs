import logging
import pathlib
import tomllib
from typing import Iterable
from urllib.parse import urlunsplit

import pandas as pd
from google.cloud import storage
from nebelung.terra_workspace import TerraWorkspace
from nebelung.utils import batch_evenly
from tqdm import tqdm

from depmap_omics_wgs.gcs import rewrite_blob

pd.set_option("display.max_columns", 30)
pd.set_option("display.max_colwidth", 50)
pd.set_option("display.max_info_columns", 30)
pd.set_option("display.max_info_rows", 20)
pd.set_option("display.max_rows", 20)
pd.set_option("display.max_seq_items", None)
pd.set_option("display.width", 200)
pd.set_option("expand_frame_repr", True)
pd.set_option("mode.chained_assignment", "warn")


def get_objects_metadata(urls: Iterable[str]) -> pd.DataFrame:
    storage_client = storage.Client(project="depmap-omics")

    blobs = {}

    with storage_client.batch(raise_exception=False):
        for url in urls:
            blob = storage.Blob.from_string(url, client=storage_client)
            bucket = storage_client.bucket(blob.bucket.name)
            blob = bucket.get_blob(blob.name)
            blobs[url] = blob

    metadata = [
        {"url": k, "crc32c": v.crc32c, "size": v.size, "gcs_obj_updated_at": v.updated}
        for k, v in blobs.items()
    ]

    return pd.DataFrame(metadata)


def get_gcs_object_metadata(urls: Iterable[str], gcp_project_id: str) -> pd.DataFrame:
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

    return df


def copy_to_cclebams(
    samples: pd.DataFrame,
    gcp_project_id: str,
    gcs_destination_bucket: str,
    gcs_destination_prefix: str,
    dry_run: bool = True,
) -> pd.DataFrame:
    """
    Copy all BAM files in the samples data frame to our bucket.

    :param samples: the data frame of samples
    :param gcp_project_id: the ID of a GCP project to use for billing
    :param gcs_destination_bucket: the name of the destination bucket
    :param gcs_destination_prefix: an object prefix for the copied BAMs
    :param dry_run: whether to skip updates to external data stores
    :return: a data frame of files we attempted to copy
    """

    print("Copying files to our bucket...")

    storage_client = storage.Client(project=gcp_project_id)

    # collect all CRAM/BAM/CRAI/BAI files to copy
    sample_files = samples.melt(
        id_vars="sample_id",
        value_vars=["internal_bai_filepath_alt", "internal_bam_filepath_alt"],
        var_name="url_kind",
        value_name="url",
    ).dropna()

    # all copied files will have same destination bucket and prefix
    dest_bucket = storage_client.bucket(
        gcs_destination_bucket, user_project=gcp_project_id
    )
    prefix = gcs_destination_prefix.strip("/")

    # keep track of copy attempts
    copy_results = []

    print(f"Checking {len(sample_files)} files to copy...")

    # can't use rewrite in a batch context, so do plain iteration
    for r in tqdm(sample_files.itertuples(index=False), total=len(sample_files)):
        url = r.url

        try:
            # construct the source blob
            src_blob = storage.Blob.from_string(url, client=storage_client)

            # construct the destination blob (named by sample ID)
            dest_file_ext = pathlib.Path(str(src_blob.name)).suffix
            dest_obj_key = "/".join([prefix, str(r.sample_id)]) + dest_file_ext
            dest_blob = storage.Blob(dest_obj_key, bucket=dest_bucket)

            # GCS rewrite operation is instantaneous if location and storage class match
            dest_blob.storage_class = src_blob.bucket.get_blob(
                src_blob.name
            ).storage_class

            new_url = urlunsplit(("gs", dest_bucket.name, dest_obj_key, "", ""))

            if dry_run:
                print(f"(skipping) Copying {dest_blob.name} to {dest_bucket.name}")
            else:
                print(f"Copying {dest_blob.name} to {dest_bucket.name}")
                rewrite_blob(src_blob, dest_blob)

            copy_results.append({"url": url, "new_url": new_url, "copied": True})

        except Exception as e:
            logging.error(f"Error copying {url} to {dest_bucket}: {e}")
            copy_results.append({"url": url, "new_url": None, "copied": False})

    return sample_files.merge(pd.DataFrame(copy_results), how="inner", on="url")


# use same config loading as when calling the module CLI
with open("config.toml", "rb") as f:
    config = tomllib.load(f)

alt_urls = pd.read_csv("~/Desktop/wes_alt.tsv", sep="\t").dropna()

workspace = TerraWorkspace(
    workspace_namespace="broad-firecloud-ccle",
    workspace_name="DepMap_WES_CN_hg38",
)

samples = workspace.get_entities("sample")[
    ["sample_id", "internal_bai_filepath", "internal_bam_filepath"]
].dropna()

samples = samples.loc[samples["sample_id"].isin(alt_urls["entity:sample_id"])]

blobs = get_objects_metadata(urls=samples["internal_bam_filepath"])

# sample_files = copy_to_cclebams(
#     samples,
#     gcp_project_id="depmap-omics",
#     gcs_destination_bucket="cclebams",
#     gcs_destination_prefix="hg38_wes",
#     dry_run=True,
# )
#
# unmelted = (
#     sample_files.pivot(index="sample_id", columns="url_kind", values="new_url")
#     .reset_index()
#     .astype("string")
#     .rename(
#         columns={
#             "internal_bai_filepath_alt": "internal_bai_filepath",
#             "internal_bam_filepath_alt": "internal_bam_filepath",
#         }
#     )
# )
#
# workspace.upload_entities(unmelted)

seq_align_updates = samples.merge(
    blobs, how="inner", left_on="internal_bam_filepath", right_on="url"
)

for _, r in seq_align_updates.iterrows():
    print(f"""
update sequencing_alignment set
    url='{r["internal_bam_filepath"]}',
    index_url='{r["internal_bai_filepath"]}',
    size={r["size"]},
    crc32c_hash='{r["crc32c"]}'
where sequencing_alignment_source='CDS' and omics_sequencing_id='{r["sample_id"]}';
""")

wes_bams = pd.read_csv("~/Desktop/wes_bams.tsv", sep="\t").dropna()
blobs = get_gcs_object_metadata(
    urls=wes_bams["internal_bam_filepath"], gcp_project_id="depmap-omics"
)
print(blobs.loc[blobs["size"].isna()])

wes_bams = wes_bams.rename(
    columns={
        "entity:sample_id": "omics_sequencing_id",
        "internal_bam_filepath": "url",
        "internal_bai_filepath": "index_url",
    }
).merge(blobs.drop(columns="created_at"), how="inner", on="url")

wes_alignments = pd.read_csv("~/Desktop/wes_alignments.tsv", sep="\t")

df = wes_bams.merge(
    wes_alignments,
    how="outer",
    on="omics_sequencing_id",
    suffixes=("", "_gumbo"),
)

df_wrong = df.loc[
    df["url"].ne(df["url_gumbo"])
    | df["index_url"].ne(df["index_url_gumbo"])
    | df["crc32c_hash"].ne(df["crc32c_hash_gumbo"])
    | df["size"].ne(df["size_gumbo"])
]

df_update = df_wrong.dropna()
df_add = df_wrong.loc[
    ~df_wrong["omics_sequencing_id"].isin(df_update["omics_sequencing_id"])
]

for _, r in df_add.iterrows():
    print(f"""
insert into sequencing_alignment (
omics_sequencing_id, index_url, url, size, crc32c_hash, reference_genome, sequencing_alignment_source
) values (
'{r["omics_sequencing_id"]}', '{r["index_url"]}', '{r["url"]}', {r["size"]}, '{r["crc32c_hash"]}', 'hg38', 'CDS'
);
""")
