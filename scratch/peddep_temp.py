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

storage_client = storage.Client(project="depmap-omics")

# collect all CRAM/BAM/CRAI/BAI files to copy
sample_files = """gs://cclebams/wgs_hg38/CDS-yb5x9U.bam
gs://cclebams/wgs_hg38/CDS-ePifjT.bam
gs://cclebams/wgs_hg38/CDS-bub9Tt.bam
gs://cclebams/wgs_hg38/CDS-pqdumK.hg38.bam
gs://cclebams/wgs_hg38/CDS-kEkc7P.hg38.bam""".splitlines()

# all copied files will have same destination bucket and prefix
dest_bucket = storage_client.bucket("peddep-share-temp", user_project="depmap-omics")

# keep track of copy attempts
copy_results = []

print(f"Checking {len(sample_files)} files to copy...")

logger = logging.getLogger()
logger.addHandler(logging.StreamHandler())
logger.setLevel(logging.INFO)

# can't use rewrite in a batch context, so do plain iteration
for url in tqdm(sample_files, total=len(sample_files)):
    try:
        # construct the source blob
        src_blob = storage.Blob.from_string(url, client=storage_client)
        dest_blob = storage.Blob(
            url.replace("cclebams", "peddep-share-temp"),
            bucket=dest_bucket,
        )

        # GCS rewrite operation is instantaneous if location and storage class match
        dest_blob.storage_class = src_blob.bucket.get_blob(src_blob.name).storage_class

        print(f"Copying {dest_blob.name} to {dest_bucket.name}")
        rewrite_blob(src_blob, dest_blob)

    except Exception as e:
        logging.error(f"Error copying {url} to {dest_bucket}: {e}")
        copy_results.append({"url": url, "new_url": None, "copied": False})
