"""
Get already-preprocessed BAM files from DepMap_WGS_CN and populate the Gumbo
`task_result` table.
"""

import os
import pathlib

import pandas as pd
from firecloud import api as firecloud_api

from gumbo_gql_client import (
    GumboClient,
    task_entity_insert_input,
    task_result_bool_exp,
    task_result_insert_input,
)
from omics_wgs_pipeline.types import CoercedDataFrame
from omics_wgs_pipeline.utils import (
    anti_join,
    df_to_model,
    expand_dict_columns,
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

    # all samples should be hg38, but this correctly assigns hg19 if necessary
    gumbo_samples["reference_assembly"] = "hg38"
    gumbo_samples["is_hg19"] = (
        gumbo_samples["cram_bam"].eq(gumbo_samples["hg19_bam_filepath"]).fillna(False)
    )
    gumbo_samples.loc[gumbo_samples["is_hg19"], "reference_assembly"] = "hg19"

    return gumbo_samples.loc[
        ~gumbo_samples["is_hg19"], ["sequencing_id", "cram_bam", "crai_bai"]
    ]


def get_wgs_cn_samples():
    wgs_cn_samples = pd.DataFrame(
        firecloud_api.get_entities(
            namespace="broad-firecloud-ccle", workspace="DepMap_WGS_CN", etype="sample"
        ).json()
    )
    wgs_cn_samples = expand_dict_columns(wgs_cn_samples)
    wgs_cn_samples = wgs_cn_samples.convert_dtypes()

    wgs_cn_samples = wgs_cn_samples.rename(
        columns={
            "name": "sequencing_id",
            "internal_bam_filepath": "bam",
            "internal_bai_filepath": "bai",
        }
    )

    wgs_cn_samples = wgs_cn_samples[["sequencing_id", "bam", "bai"]].dropna()

    return wgs_cn_samples


def populate_task_entities(task_entities, task_results):
    new_task_entities = task_results.loc[
        ~task_results["sequencing_id"].isin(task_entities["sequencing_id"]),
        ["sequencing_id"],
    ]

    if len(new_task_entities) > 0:
        gumbo_client.insert_task_entities(
            username="dmccabe",
            objects=df_to_model(new_task_entities, task_entity_insert_input),
        )

        task_entities = model_to_df(
            gumbo_client.sequencing_task_entities(),
            CoercedDataFrame,
            remove_unknown_cols=False,
        )

    return task_entities


def make_task_results(task_results):
    task_results = task_results.merge(
        task_entities.rename(columns={"id": "task_entity_id"}),
        how="left",
        on="sequencing_id",
    )
    task_results = task_results.loc[:, ["task_entity_id", "bam", "bai"]]

    task_results["workflow_name"] = "preprocess_wgs_sample"
    task_results["task_name"] = "gatk_applybqsr"
    task_results["created_at"] = pd.Timestamp.now(tz="UTC").isoformat()

    task_results = task_results.melt(
        id_vars=["task_entity_id", "workflow_name", "task_name", "created_at"],
        value_vars=["bam", "bai"],
        var_name="label",
        value_name="url",
    )

    task_results["format"] = task_results["url"].apply(
        lambda x: pathlib.Path(x).suffix[1:].upper()
    )

    return task_results


gumbo_samples = get_gumbo_samples()
wgs_cn_samples = get_wgs_cn_samples()

task_results = wgs_cn_samples.loc[
    wgs_cn_samples["sequencing_id"].isin(gumbo_samples["sequencing_id"])
]

task_entities = model_to_df(
    gumbo_client.sequencing_task_entities(), CoercedDataFrame, remove_unknown_cols=False
)
task_entities = populate_task_entities(task_entities, task_results)

task_results = make_task_results(task_results)

existing_task_results = model_to_df(
    gumbo_client.get_task_results(
        task_result_bool_exp(
            workflow_name={"eq": "preprocess_wgs_sample"},
            task_name={"eq": "gatk_applybqsr"},
        )
    ),
    CoercedDataFrame,
    remove_unknown_cols=False,
)

existing_task_results = expand_dict_columns(
    existing_task_results, name_columns_with_parent=False
)

new_task_results = anti_join(
    task_results,
    existing_task_results,
    on=["workflow_name", "task_name", "label", "url", "format"],
)

if len(new_task_results) > 0:
    gumbo_client.insert_task_results(
        username="dmccabe",
        objects=df_to_model(new_task_results, task_result_insert_input),
    )
