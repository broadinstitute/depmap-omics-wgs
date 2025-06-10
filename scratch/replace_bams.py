import logging
import pathlib
from urllib.parse import urlunsplit

import pandas as pd
import tqdm
from google.cloud import storage
from nebelung.terra_workspace import TerraWorkspace

from depmap_omics_wgs.gcs import copy_to_cclebams, get_objects_metadata

pd.set_option("display.max_columns", 30)
pd.set_option("display.max_colwidth", 50)
pd.set_option("display.max_info_columns", 30)
pd.set_option("display.max_info_rows", 20)
pd.set_option("display.max_rows", 20)
pd.set_option("display.max_seq_items", None)
pd.set_option("display.width", 200)
pd.set_option("expand_frame_repr", True)
pd.set_option("mode.chained_assignment", "warn")

logger = logging.getLogger()
logger.addHandler(logging.StreamHandler())
logger.setLevel(logging.INFO)

ws = TerraWorkspace("broad-firecloud-ccle", "depmap-omics-wgs")
samples = ws.get_entities("sample")
sample_sets = ws.get_entities("sample_set")
sample_ids = [
    x["entityName"]
    for x in sample_sets.loc[sample_sets["sample_set_id"].eq("25q2"), "samples"].values[
        0
    ]["items"]
]
samples = samples.loc[
    samples["sample_id"].isin(sample_ids),
    ["sample_id", "analysis_ready_bam", "analysis_ready_bai"],
]

copied_files = copy_to_cclebams(
    samples,
    gcp_project_id="depmap-omics",
    gcs_destination_bucket="cclebams",
    gcs_destination_prefix="wgs_hg38",
    overwrite=True,
    dry_run=False,
)

om = get_objects_metadata(copied_files["new_url"])
om = om.loc[om["url"].str.endswith(".bam")]
om.to_csv("~/Desktop/om.csv", index=False)

samples_updated = samples.copy()

for c in ["analysis_ready_bai", "analysis_ready_bam"]:
    sample_file_urls = copied_files.loc[copied_files["copied"], ["url", "new_url"]]
    samples_updated = samples_updated.merge(
        sample_file_urls, how="left", left_on=c, right_on="url"
    )
    samples_updated[c] = samples_updated["new_url"]
    samples_updated = samples_updated.drop(columns=["url", "new_url"])

ws.upload_entities(samples_updated)

legacy_terra_workspace = TerraWorkspace("broad-firecloud-ccle", "DepMap_WGS_CN")
legacy_terra_workspace.upload_entities(
    samples_updated.rename(
        columns={
            "analysis_ready_bam": "internal_bam_filepath",
            "analysis_ready_bai": "internal_bai_filepath",
        }
    )
)
