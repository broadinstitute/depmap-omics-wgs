import logging
import pathlib
from typing import Iterable, Optional
from urllib.parse import urlunsplit

import pandas as pd
from google.cloud import storage
from nebelung.utils import type_data_frame
from pandera.typing import DataFrame as TypedDataFrame
from tqdm.asyncio import tqdm

from depmap_omics_wgs.types import (
    AlignedSamplesWithObjectMetadata,
    CopiedSampleFiles,
    GcsObject,
)


def get_objects_metadata(urls: Iterable[str]) -> TypedDataFrame[GcsObject]:
    """
    Get metadata (size, CRC32 hash, and updated_at timestamp) for a list of GCS URLs.

    :param urls: iterable of GCS URLs
    :return: data frame of metadata
    """

    storage_client = storage.Client()

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

    df = pd.DataFrame(metadata)

    # batching without raising exceptions makes all columns NA if a file was missing
    df = df.dropna()
    df = df.astype({"gcs_obj_updated_at": "datetime64[ns, UTC]"})

    # drop time component
    df["gcs_obj_updated_at"] = df["gcs_obj_updated_at"].dt.date.astype("str")

    return type_data_frame(df, GcsObject)


def copy_to_cclebams(
    samples: TypedDataFrame[AlignedSamplesWithObjectMetadata],
    gcp_project_id: str,
    gcs_destination_bucket: str,
    gcs_destination_prefix: str,
    overwrite: bool = False,
    dry_run: bool = True,
) -> TypedDataFrame[CopiedSampleFiles]:
    """
    Copy all BAM files in the samples data frame to our bucket.

    :param samples: the data frame of samples
    :param gcp_project_id: the ID of a GCP project to use for billing
    :param gcs_destination_bucket: the name of the destination bucket
    :param gcs_destination_prefix: an object prefix for the copied BAMs
    :param overwrite: whether to overwrite existing files at the destination
    :param dry_run: whether to skip updates to external data stores
    :return: a data frame of files we attempted to copy
    """

    logging.info("Copying files to our bucket...")

    storage_client = storage.Client(project=gcp_project_id)

    # collect all CRAM/BAM/CRAI/BAI files to copy
    sample_files = samples.melt(
        id_vars="sample_id",
        value_vars=["analysis_ready_bai", "analysis_ready_bam"],
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

    logging.info(f"Checking {len(sample_files)} files to copy...")

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

            if dest_blob.exists() and not overwrite:
                logging.info(f"{dest_blob.name} already exists in {dest_bucket.name}")
                copy_results.append({"url": url, "new_url": new_url, "copied": True})

            else:
                if dry_run:
                    logging.info(
                        f"(skipping) Copying {dest_blob.name} to {dest_bucket.name}"
                    )
                else:
                    logging.info(f"Copying {dest_blob.name} to {dest_bucket.name}")
                    rewrite_blob(src_blob, dest_blob)

                copy_results.append({"url": url, "new_url": new_url, "copied": True})

        except Exception as e:
            logging.error(f"Error copying {url} to {dest_bucket}: {e}")
            copy_results.append({"url": url, "new_url": None, "copied": False})

    sample_files = sample_files.merge(pd.DataFrame(copy_results), how="inner", on="url")

    return type_data_frame(sample_files, CopiedSampleFiles)


def rewrite_blob(src_blob: storage.Blob, dest_blob: storage.Blob) -> None:
    """
    Use Blob.rewrite to copy blob to another location.

    :param src_blob: a storage Blob to copy from
    :param dest_blob: a storage Blob to copy to
    """

    logging.info(f"Copying {dest_blob.name}...")

    # assume the rewrite will take multiple requests and keep checking as long as
    # the request returns a rewrite token
    token: Optional[str] = None

    while True:
        token, bytes_rewritten, total_bytes = dest_blob.rewrite(
            source=src_blob, token=token, timeout=60
        )

        logging.info(
            "Copied {bytes_rewritten} / {total_bytes} bytes ({prop}%)".format(
                bytes_rewritten=bytes_rewritten,
                total_bytes=total_bytes,
                prop=round(bytes_rewritten / total_bytes * 100),
            )
        )

        if token is None:
            break
