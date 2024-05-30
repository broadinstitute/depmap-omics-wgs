import os

import numpy as np
import pandas as pd
from firecloud import api as firecloud_api
from google.cloud import storage

from omics_wgs_pipeline.utils import expand_dict_columns

pd.set_option("display.max_columns", 30)
pd.set_option("display.max_colwidth", 50)
pd.set_option("display.max_info_columns", 30)
pd.set_option("display.max_info_rows", 20)
pd.set_option("display.max_rows", 20)
pd.set_option("display.max_seq_items", None)
pd.set_option("display.width", 200)
pd.set_option("expand_frame_repr", True)
pd.set_option("mode.chained_assignment", "warn")

# get the absolute path to the root directory of the repo
root_dir = os.path.abspath(os.path.join(os.path.dirname(__name__)))


def get_terra_outputs(namespace, workspace):
    wgs_cn_samples = pd.DataFrame(
        firecloud_api.get_entities(
            namespace=namespace, workspace=workspace, etype="sample"
        ).json()
    )
    wgs_cn_samples = expand_dict_columns(wgs_cn_samples, name_columns_with_parent=False)
    wgs_cn_samples = wgs_cn_samples.convert_dtypes()
    wgs_cn_samples = wgs_cn_samples.rename(columns={"name": "sample_id"})

    return (
        wgs_cn_samples[
            [
                "sample_id",
                "msisensor2_score",
                "msisensor2_output",
                "msisensor2_output_dis",
                "msisensor2_output_somatic",
            ]
        ]
        .dropna()
        .astype(
            {
                "sample_id": "string",
                "msisensor2_score": "float64",
                "msisensor2_output": "string",
                "msisensor2_output_dis": "string",
                "msisensor2_output_somatic": "string",
            }
        )
    )


def save_msi_outputs(storage_client, df, old_or_new):
    local_dir = os.path.join(root_dir, "data", old_or_new)

    for _, r in df.iterrows():
        for c in [
            "msisensor2_output",
            "msisensor2_output_dis",
            "msisensor2_output_somatic",
        ]:
            sample_dir = os.path.join(local_dir, r["sample_id"])
            os.makedirs(sample_dir, exist_ok=True)
            local_file = os.path.join(sample_dir, ".".join([r["sample_id"], c]))

            blob = storage.Blob.from_string(r[c], client=storage_client)
            blob.download_to_filename(local_file)

    pass


def compare_scores(old_samples, new_samples):
    score_comp = old_samples[["sample_id", "msisensor2_score"]].merge(
        new_samples[["sample_id", "msisensor2_score"]],
        on="sample_id",
        how="inner",
        suffixes=("_old", "_new"),
    )
    score_med_abs_dev = np.abs(
        score_comp["msisensor2_score_new"] - score_comp["msisensor2_score_old"]
    ).median()
    score_max_abs_dev = np.abs(
        score_comp["msisensor2_score_new"] - score_comp["msisensor2_score_old"]
    ).max()
    print(f"Median absolute deviation: {round(score_med_abs_dev, 2)}")
    print(f"Max absolute deviation: {round(score_max_abs_dev, 2)}")


old_samples = get_terra_outputs(
    namespace="broad-firecloud-ccle", workspace="DepMap_WGS_CN"
)
new_samples = get_terra_outputs(
    namespace="broad-firecloud-ccle", workspace="omics_wgs_pipeline"
)

sample_ids = set(old_samples["sample_id"]).intersection(set(new_samples["sample_id"]))
old_samples = old_samples.loc[old_samples["sample_id"].isin(sample_ids)]
new_samples = new_samples.loc[new_samples["sample_id"].isin(sample_ids)]

storage_client = storage.Client(project="depmap-omics")
save_msi_outputs(storage_client, old_samples, "old")
save_msi_outputs(storage_client, new_samples, "new")

compare_scores(old_samples, new_samples)
