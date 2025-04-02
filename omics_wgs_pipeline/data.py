import datetime
import json
import logging

import pandas as pd
import requests
from firecloud import api as firecloud_api
from nebelung.terra_workspace import TerraWorkspace
from nebelung.utils import call_firecloud_api, type_data_frame

from gumbo_gql_client import task_entity_insert_input, task_result_insert_input
from omics_wgs_pipeline.types import (
    GumboClient,
    GumboTaskEntity,
    GumboTaskResult,
    GumboWgsSequencing,
    TerraSample,
    TypedDataFrame,
)
from omics_wgs_pipeline.utils import (
    compute_uuidv3,
    get_gcs_object_metadata,
    model_to_df,
)


def do_refresh_terra_samples(
    terra_workspace: TerraWorkspace,
    gumbo_client: GumboClient,
    ref_urls: dict[str, dict[str, str]],
) -> None:
    """
    Update the Terra `sample` data table using ground truth data from Gumbo.

    :param terra_workspace: a `TerraWorkspace` instance
    :param gumbo_client: a `GumboClient` instance
    :param ref_urls: a nested dictionary of genomes and their reference file URLs
    """

    # get long data frame of both GP-delivered and CDS (analysis ready) CRAM/BAMs
    wgs_sequencings = model_to_df(
        gumbo_client.wgs_sequencing_alignments(), GumboWgsSequencing
    )

    # make wide, separating delivery and analysis-ready CRAM/BAMs
    samples = (
        wgs_sequencings.loc[wgs_sequencings["sequencing_alignment_source"].eq("GP")]
        .drop(columns="sequencing_alignment_source")
        .rename(
            columns={
                "omics_sequencing_id": "sample_id",
                "sequencing_alignment_id": "delivery_sequencing_alignment_id",
                "url": "delivery_cram_bam",
                "index_url": "delivery_crai_bai",
                "reference_genome": "delivery_ref",
            }
        )
        .merge(
            wgs_sequencings.loc[
                wgs_sequencings["sequencing_alignment_source"].eq("CDS")
            ]
            .drop(columns="sequencing_alignment_source")
            .rename(
                columns={
                    "omics_sequencing_id": "sample_id",
                    "sequencing_alignment_id": "aligned_sequencing_alignment_id",
                    "url": "analysis_ready_bam",
                    "index_url": "analysis_ready_bai",
                    "reference_genome": "ref",
                }
            ),
            how="outer",
            on="sample_id",
        )
    )

    samples["delivery_file_format"] = (
        samples["delivery_cram_bam"].str.rsplit(".", n=1).str.get(1).str.upper()
    )

    # set reference genome columns
    samples = set_ref_urls(samples, ref_urls)

    # validate types
    samples = type_data_frame(samples, TerraSample)

    # delete obsolete samples (e.g. ones that have been blacklisted since the last sync)
    terra_samples = terra_workspace.get_entities("sample")
    terra_workspace.delete_entities(
        entity_type="sample",
        entity_ids=set(terra_samples["sample_id"]).difference(
            set(samples["sample_id"])
        ),
    )

    sample_ids = samples.pop("sample_id")
    samples.insert(0, "entity:sample_id", sample_ids)
    terra_workspace.upload_entities(df=samples)


def set_ref_urls(
    samples: pd.DataFrame, ref_urls: dict[str, dict[str, str]]
) -> pd.DataFrame:
    """
    Populate columns in the `samples` data frame with URLs for reference genome files.

    :param samples: the data frame of sample data
    :param ref_urls: a dictionary of reference genome names (e.g. "hg38") and URLs of
    files
    :return: the `samples` data frame with columns for reference genome URLs
    """

    # set default reference for (re)alignment
    samples["ref"] = samples["ref"].fillna("hg38")

    # join reference URLs to for `delivery_` and `analysis_ready_` CRAM/BAMs
    ref_df = pd.DataFrame(ref_urls.values())
    ref_df["ref"] = ref_urls.keys()
    samples = samples.merge(ref_df, how="left", on="ref")

    ref_df.columns = "delivery_" + ref_df.columns

    return samples.merge(ref_df, how="left", on="delivery_ref")


def do_refresh_legacy_terra_samples(
    terra_workspace: TerraWorkspace,
    legacy_terra_workspace: TerraWorkspace,
    sample_set_id: str,
    gumbo_client: GumboClient,
) -> None:
    """
    Update the sample data table in the legacy production Terra workspace.

    :param terra_workspace: the Terra workspace where new workflows have run
    :param legacy_terra_workspace: the legacy Terra workspace where new workflow outputs
    need to be synced to
    :param sample_set_id: the ID for a sample set to create or update for unprocessed
    samples
    :param gumbo_client: a `GumboClient` instance
    """

    # identify samples in legacy workspace that have new analysis-ready BAMs
    samples = terra_workspace.get_entities("sample")
    legacy_samples = legacy_terra_workspace.get_entities("sample")

    src_samples = samples[
        [
            "sample_id",
            "analysis_ready_bam",
            "analysis_ready_bai",
            "cnv_cn_by_gene_weighted_mean",
            "cnv_segments",
        ]
    ].rename(
        columns={
            "analysis_ready_bam": "internal_bam_filepath",
            "analysis_ready_bai": "internal_bai_filepath",
        }
    )

    if "internal_bam_filepath" not in legacy_samples.columns:
        legacy_samples["internal_bam_filepath"] = pd.NA

    if "internal_bai_filepath" not in legacy_samples.columns:
        legacy_samples["internal_bai_filepath"] = pd.NA

    if "cnv_cn_by_gene_weighted_mean" not in legacy_samples.columns:
        legacy_samples["cnv_cn_by_gene_weighted_mean"] = pd.NA

    if "cnv_segments" not in legacy_samples.columns:
        legacy_samples["cnv_segments"] = pd.NA

    dest_samples = legacy_samples[
        [
            "sample_id",
            "internal_bam_filepath",
            "internal_bai_filepath",
            "cnv_cn_by_gene_weighted_mean",
            "cnv_segments",
        ]
    ]

    diff = src_samples.merge(
        dest_samples, on="sample_id", how="left", suffixes=("_src", "_dest")
    )

    to_upsert = diff.loc[
        (
            diff["internal_bam_filepath_src"].notna()
            & diff["internal_bam_filepath_dest"].isna()
        )
        | (
            diff["internal_bai_filepath_src"].notna()
            & diff["internal_bai_filepath_dest"].isna()
        )
        | (
            diff["cnv_cn_by_gene_weighted_mean_src"].notna()
            & diff["cnv_cn_by_gene_weighted_mean_dest"].isna()
        )
        | (diff["cnv_segments_src"].notna() & diff["cnv_segments_dest"].isna()),
        [
            "sample_id",
            "internal_bam_filepath_src",
            "internal_bai_filepath_src",
            "cnv_cn_by_gene_weighted_mean_src",
            "cnv_segments_src",
        ],
    ].rename(
        columns={
            "internal_bam_filepath_src": "internal_bam_filepath",
            "internal_bai_filepath_src": "internal_bai_filepath",
            "cnv_cn_by_gene_weighted_mean_src": "cnv_cn_by_gene_weighted_mean",
            "cnv_segments_src": "cnv_segments",
        }
    )

    if len(to_upsert) > 0:
        legacy_terra_workspace.upload_entities(to_upsert)
        legacy_samples = legacy_terra_workspace.get_entities("sample")

    # get list of unprocessesd sample/sequencing IDs and create/update a sample set
    res = gumbo_client.unprocessed_sequencings()

    unprocessed_sample_ids = [
        x.omics_sequencing_id
        for x in res.records
        if x.omics_sequencing_id in legacy_samples["sample_id"].values
    ]

    try:
        # update existing sample set
        sample_set = call_firecloud_api(
            firecloud_api.get_entity,
            namespace=legacy_terra_workspace.workspace_namespace,
            workspace=legacy_terra_workspace.workspace_name,
            etype="sample_set",
            ename=sample_set_id,
        )

        sample_set["attributes"]["samples"]["items"] = [
            {"entityType": "sample", "entityName": x} for x in unprocessed_sample_ids
        ]

        _ = call_firecloud_api(
            firecloud_api.update_entity,
            namespace=legacy_terra_workspace.workspace_namespace,
            workspace=legacy_terra_workspace.workspace_name,
            etype="sample_set",
            ename=sample_set_id,
            updates=[
                {
                    "op": "AddUpdateAttribute",
                    "attributeName": "samples",
                    "addUpdateAttribute": sample_set["attributes"]["samples"],
                }
            ],
        )

    except requests.exceptions.RequestException as e:
        if "404" not in str(e):
            raise e

        # create new sample set
        legacy_terra_workspace.create_entity_set(
            entity_type="sample",
            entity_ids=unprocessed_sample_ids,
            entity_set_id=sample_set_id,
        )


def pick_best_task_results(
    task_results: TypedDataFrame[GumboTaskResult],
) -> pd.DataFrame:
    """
    Pick the most recent task result (files and values) for each sample+label
    combination and make a wide data frame for use in a Terra `sample` data table.

    :param task_results: a data frame of Gumbo task results
    :return: a wide data frame of samples and their outputs
    """

    # check that distinct output labels (i.e. data table columns) don't come from
    # different workflows
    workflow_names_per_label = (
        task_results[["workflow_name", "label"]]
        .drop_duplicates()["label"]
        .value_counts()
    )

    bad_labels = workflow_names_per_label.loc[workflow_names_per_label.gt(1)]

    if len(bad_labels) > 0:
        raise ValueError(
            "\n".join(
                [
                    "Cannot sync task results back to Terra. "
                    "Inconsistent workflow name and label combinations:",
                    str(bad_labels.to_frame().reset_index()),
                ]
            )
        )

    # get unique file (URL) outputs for each sample+label
    task_result_files = (
        task_results[["sample_id", "label", "url", "workflow_version", "completed_at"]]
        .dropna(subset=["sample_id", "label", "url"])
        .drop_duplicates(subset=["sample_id", "label", "url"])
    )

    # pick the most recent file based on workflow version and date
    best_task_result_files = (
        task_result_files.sort_values(
            ["workflow_version", "completed_at"], ascending=False
        )
        .groupby(["sample_id", "label"])
        .nth(0)
        .pivot(index="sample_id", columns="label", values="url")
        .reset_index()
    )

    # get unique value outputs for each sample+label
    task_result_values = task_results[
        ["sample_id", "label", "value", "workflow_version", "completed_at"]
    ].dropna(subset=["sample_id", "label", "value"])

    # drop NA values
    task_result_values = task_result_values.loc[
        task_result_values["value"].ne({"value": None})
    ]

    # need to temporarily convert values to JSON strings so that they can be de-deduped
    task_result_values["value"] = task_result_values["value"].apply(json.dumps)

    task_result_values = task_result_values.drop_duplicates(
        subset=["sample_id", "label", "value"]
    )

    task_result_values["value"] = task_result_values["value"].apply(json.loads)

    # pick the most recent value based on workflow version and date
    best_task_result_values = (
        task_result_values.sort_values(
            ["workflow_version", "completed_at"], ascending=False
        )
        .groupby(["sample_id", "label"])
        .nth(0)
        .pivot(index="sample_id", columns="label", values="value")
        .reset_index()
        .set_index("sample_id")
        .fillna(pd.NA)
    )

    # extract the values from the dictionaries
    best_task_result_values = (
        best_task_result_values.map(lambda x: x["value"] if x is not pd.NA else pd.NA)
        .reset_index()
        .fillna(pd.NA)
    )

    # join files and values and guess their dtypes
    best_task_results = best_task_result_files.merge(
        best_task_result_values, how="outer", on="sample_id"
    ).convert_dtypes()

    return best_task_results


def put_task_results(
    gumbo_client: GumboClient,
    terra_workspace: TerraWorkspace,
    gcp_project_id: str,
    uuid_namespace: str,
    since: datetime.datetime | None = None,
) -> None:
    """
    Upsert a list of new `task_result` records to Gumbo and associate them with a new
    `terra_sync` record.

    :param terra_workspace: a `TerraWorkspace` instance
    :param gumbo_client: a `GumboClient` instance
    :param gcp_project_id: the ID of a GCP project to use for billing
    :param uuid_namespace: a namespace for generated a UUIDv3 ID for each record
    :param since: don't collect outputs for job submissions before this `datetime`
    """

    # get all current omics_sequencing IDs
    valid_seq_ids = set(x.sequencing_id for x in gumbo_client.sequencing_ids().records)

    # get the outputs from the workspace's completed jobs
    all_outputs = terra_workspace.collect_workflow_outputs(since)
    outputs = []

    for x in all_outputs:
        if (
            x.terra_entity_type != "sample"
            or "sample_id" not in x.terra_workflow_inputs
            or x.terra_workflow_inputs["sample_id"] not in valid_seq_ids
        ):
            # not persisting outputs for jobs that didn't operate on samples or for
            # unknown sequencing/sample IDs
            continue

        o = x.model_copy()

        # replace the entity name with the more reliable `sample_id` workflow input
        o.terra_entity_name = o.terra_workflow_inputs["sample_id"]

        # do trivial conversion from `TaskResult` to `task_result_insert_input` type
        outputs.append(task_result_insert_input(**o.model_dump()))

    logging.info("Getting existing task entity records for sequencings")
    task_entities = model_to_df(
        gumbo_client.sequencing_task_entities(),
        GumboTaskEntity,
        remove_unknown_cols=True,
    )

    # check if any need to be created (a new `task_result` record must belong to one)
    req_seq_ids = set([x.terra_entity_name for x in outputs])
    missing_seq_ids = req_seq_ids.intersection(valid_seq_ids).difference(
        task_entities["sequencing_id"]
    )

    if len(missing_seq_ids) > 0:
        logging.info(f"Creating {len(missing_seq_ids)} missing task entity records")
        gumbo_client.insert_task_entities(
            username=gumbo_client.username,
            objects=[
                task_entity_insert_input(sequencing_id=x) for x in missing_seq_ids
            ],
        )

        # get the records again
        task_entities = model_to_df(
            gumbo_client.sequencing_task_entities(), GumboTaskEntity
        )

    # create a mapping from sequencing to task entity IDs
    task_entities = task_entities.set_index("sequencing_id").rename(
        columns={"id": "task_entity_id"}
    )

    # get GCS object metadata for this set of output URLs
    object_metadata = get_gcs_object_metadata(
        [str(x.url) for x in outputs if x.url is not None], gcp_project_id
    ).set_index("url")

    for i in range(len(outputs)):
        # assign final missing task result values
        outputs[i].task_entity_id = task_entities.loc[
            outputs[i].terra_entity_name, "task_entity_id"
        ]

        if outputs[i].url in object_metadata.index:
            om = object_metadata.loc[outputs[i].url]
            outputs[i].size = int(om["size"])
            outputs[i].crc_32_c_hash = om["crc32c_hash"]

        # assign a persistent UUID based on the record's values
        outputs[i].id = compute_uuidv3(
            outputs[i].model_dump(mode="json", by_alias=True),
            uuid_namespace,
            # use all fields with known values
            keys={
                "crc32c_hash",
                "completed_at",
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

    logging.info(f"Upserting {len(outputs)} task results into Gumbo")
    res = gumbo_client.insert_terra_sync(
        username=gumbo_client.username,
        terra_workspace_namespace=str(outputs[0].terra_workspace_namespace),
        terra_workspace_name=str(outputs[0].terra_workspace_name),
        task_results=outputs,
    )

    sync_id = res.insert_terra_sync.returning[0].id  # pyright: ignore
    logging.info(f"Created Terra sync record {sync_id}")
