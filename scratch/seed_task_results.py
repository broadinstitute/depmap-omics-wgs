"""
Get already-preprocessed outputs from DepMap_WGS_CN and populate the Gumbo `task_result`
table.
"""

import os
import pathlib

import pandas as pd
from firecloud import api as firecloud_api

from gumbo_gql_client import task_entity_insert_input, task_result_insert_input
from omics_wgs_pipeline.types import CoercedDataFrame, GumboClient, GumboTaskEntity
from omics_wgs_pipeline.utils import (
    batch_evenly,
    compute_uuidv3,
    df_to_model,
    expand_dict_columns,
    get_gcs_object_metadata,
    model_to_df,
)

pd.set_option("display.max_columns", 30)
pd.set_option("display.max_colwidth", 50)
pd.set_option("display.max_info_columns", 30)
pd.set_option("display.max_info_rows", 20)
pd.set_option("display.max_rows", 20)
pd.set_option("display.max_seq_items", None)
pd.set_option("display.width", 200)
pd.set_option("expand_frame_repr", True)
pd.set_option("mode.chained_assignment", "warn")

gumbo_client = GumboClient(
    url=os.environ["HASURA_URL"],
    headers={"X-Hasura-Admin-Secret": os.environ["HASURA_ADMIN_SECRET"]},
    username="dmccabe",
)


def get_gumbo_samples():
    gumbo_samples = model_to_df(
        gumbo_client.wgs_sequencings(),
        CoercedDataFrame,
        remove_unknown_cols=False,
        mutator=lambda df: df.rename(
            columns={
                "hg_19_bai_filepath": "hg19_bai_filepath",
                "hg_19_bam_filepath": "hg19_bam_filepath",
                "hg_38_crai_filepath": "hg38_crai_filepath",
                "hg_38_cram_filepath": "hg38_cram_filepath",
            }
        ),
    )

    gumbo_samples = gumbo_samples.convert_dtypes()

    gumbo_samples["cram_bam"] = (
        gumbo_samples["hg38_cram_filepath"]
        .fillna(gumbo_samples["bam_filepath"])
        .fillna(gumbo_samples["hg19_bam_filepath"])
    )

    gumbo_samples["crai_bai"] = (
        gumbo_samples["hg38_crai_filepath"]
        .fillna(gumbo_samples["bai_filepath"])
        .fillna(gumbo_samples["hg19_bai_filepath"])
    )

    # all new samples should be hg38, but this removes hg19 ones if necessary
    gumbo_samples["reference_assembly"] = "hg38"
    gumbo_samples["is_hg19"] = (
        gumbo_samples["cram_bam"].eq(gumbo_samples["hg19_bam_filepath"]).fillna(False)
    )
    gumbo_samples = gumbo_samples.loc[~gumbo_samples["is_hg19"]]

    return gumbo_samples.loc[
        ~gumbo_samples["is_hg19"], ["sequencing_id", "cram_bam", "crai_bai"]
    ]


def get_wgs_cn_samples():
    wgs_cn_samples = pd.DataFrame(
        firecloud_api.get_entities(
            namespace="broad-firecloud-ccle", workspace="DepMap_WGS_CN", etype="sample"
        ).json()
    )
    wgs_cn_samples = expand_dict_columns(wgs_cn_samples, name_columns_with_parent=False)
    wgs_cn_samples = wgs_cn_samples.convert_dtypes()

    wgs_cn_samples = wgs_cn_samples.rename(
        columns={
            "name": "sequencing_id",
            "internal_bam_filepath": "bam",
            "internal_bai_filepath": "bai",
        }
    )

    return wgs_cn_samples


def populate_task_entities(task_entities, task_results):
    new_task_entities = task_results.loc[
        ~task_results["sequencing_id"].isin(task_entities["sequencing_id"]),
        ["sequencing_id"],
    ]

    if len(new_task_entities) == 0:
        return task_entities

    gumbo_client.insert_task_entities(
        username="dmccabe",
        objects=df_to_model(new_task_entities, task_entity_insert_input),
    )

    return model_to_df(gumbo_client.sequencing_task_entities(), GumboTaskEntity)


def make_task_results(task_results):
    task_results = task_results.merge(
        task_entities.rename(columns={"id": "task_entity_id"}),
        how="left",
        on="sequencing_id",
    )

    task_results = task_results.loc[
        :,
        [
            "task_entity_id",
            "bam",
            "bai",
            "msisensor2_score",
            "msisensor2_output",
            "msisensor2_output_dis",
            "msisensor2_output_somatic",
        ],
    ]

    # preprocess_wgs_sample
    task_results_prep = task_results.copy(deep=True)
    task_results_prep = task_results_prep.loc[
        :, ["task_entity_id", "bam", "bai"]
    ].dropna()
    task_results_prep["workflow_name"] = "preprocess_wgs_sample"

    task_results_prep = task_results_prep.melt(
        id_vars=["task_entity_id", "workflow_name"],
        value_vars=["bam", "bai"],
        var_name="label",
        value_name="url",
    )

    task_results_prep["format"] = task_results_prep["url"].apply(
        lambda x: pathlib.Path(x).suffix[1:].upper()
    )

    # infer_msi_status
    task_results_msi = task_results.copy(deep=True)
    task_results_msi = task_results_msi.loc[
        :,
        [
            "task_entity_id",
            "msisensor2_score",
            "msisensor2_output",
            "msisensor2_output_dis",
            "msisensor2_output_somatic",
        ],
    ].dropna()
    task_results_msi["workflow_name"] = "infer_msi_status"

    task_results_msi_urls = task_results_msi.melt(
        id_vars=["task_entity_id", "workflow_name"],
        value_vars=[
            "msisensor2_output",
            "msisensor2_output_dis",
            "msisensor2_output_somatic",
        ],
        var_name="label",
        value_name="url",
    )

    task_results_msi_values = task_results_msi.melt(
        id_vars=["task_entity_id", "workflow_name"],
        value_vars=["msisensor2_score"],
        var_name="label",
        value_name="value",
    )
    task_results_msi_values["value"] = task_results_msi_values["value"].apply(
        lambda x: {"value": x}
    )

    all_task_results = pd.concat(
        [task_results_prep, task_results_msi_urls, task_results_msi_values]
    )
    all_task_results.loc[:, ["url", "format", "value"]] = all_task_results.loc[
        :, ["url", "format", "value"]
    ].fillna(pd.NA)

    return all_task_results


def make_task_result_input(x):
    d = x.loc[~x.isna()].to_dict()
    o = task_result_insert_input(**d)

    o.id = compute_uuidv3(
        o.model_dump(mode="json", by_alias=True),
        "00000000-0000-0000-0000-000000000000",
        # use all fields with known values
        keys={
            "crc32c_hash",
            "created_at",
            "format",
            "label",
            "size",
            "task_entity_id",
            "terra_entity_name",
            "terra_entity_type",
            "terra_method_config_name",
            "terra_method_config_namespace",
            "terra_submission_id",
            "terra_workflow_id",
            "terra_workflow_inputs",
            "terra_workflow_root_dir",
            "terra_workspace_id",
            "terra_workspace_name",
            "terra_workspace_namespace",
            "url",
            "value",
            "workflow_name",
            "workflow_source_url",
            "workflow_version",
        },
    )

    return o


gumbo_samples = get_gumbo_samples()
wgs_cn_samples = get_wgs_cn_samples()

task_results = wgs_cn_samples.loc[
    wgs_cn_samples["sequencing_id"].isin(gumbo_samples["sequencing_id"])
]

task_entities = model_to_df(gumbo_client.sequencing_task_entities(), GumboTaskEntity)
task_entities = populate_task_entities(task_entities, task_results)

task_results = make_task_results(task_results)

object_metadata = get_gcs_object_metadata(task_results["url"].dropna(), "depmap-omics")
task_results = task_results.merge(object_metadata, how="left", on="url")

# set `created_at` to beginnning of UNIX epoch for values since we don't have them from
# GCS
task_results["created_at"] = task_results["created_at"].fillna(pd.Timestamp(0))

task_result_inserts = [make_task_result_input(x) for _, x in task_results.iterrows()]

for batch in batch_evenly(task_result_inserts, 500):
    gumbo_client.insert_task_results(username=gumbo_client.username, objects=batch)
