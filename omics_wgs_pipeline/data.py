import datetime
import json
import pathlib

import pandas as pd
from click import echo
from nebelung.terra_workflow import TerraWorkflow
from nebelung.terra_workspace import TerraWorkspace
from nebelung.utils import expand_dict_columns, type_data_frame

from gumbo_gql_client import (
    task_entity_insert_input,
    task_result_bool_exp,
    task_result_insert_input,
)
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


def make_terra_samples(
    gumbo_client: GumboClient, ref_urls: dict[str, dict[str, str]]
) -> TypedDataFrame[TerraSample]:
    """
    Make a data frame to use as a Terra `sample` data table using ground truth data from
    Gumbo.

    :param gumbo_client: a `GumboClient` instance
    :param ref_urls: a nested dictionary of genomes and their reference file URLs
    :return: a data frame to use as a Terra `sample` data table
    """

    wgs_sequencings = model_to_df(
        gumbo_client.wgs_sequencings(),
        GumboWgsSequencing,
        remove_unknown_cols=True,
        mutator=lambda df: df.rename(
            columns={
                "hg_19_bai_filepath": "hg19_bai_filepath",
                "hg_19_bam_filepath": "hg19_bam_filepath",
                "hg_38_crai_filepath": "hg38_crai_filepath",
                "hg_38_cram_filepath": "hg38_cram_filepath",
            }
        ),
    )

    task_results = model_to_df(
        gumbo_client.get_task_results(
            task_result_bool_exp(
                workflow_name={
                    "in_": ["preprocess_wgs_sample", "infer_msi_status"],
                },  # pyright: ignore
            )
        ),
        GumboTaskResult,
        remove_unknown_cols=True,
        mutator=lambda df: expand_dict_columns(df, except_cols=["value"]).rename(
            columns={"task_entity__sequencing_id": "sample_id", "value__value": "value"}
        ),
    )

    return join_existing_results_to_samples(wgs_sequencings, task_results, ref_urls)


def join_existing_results_to_samples(
    wgs_sequencings: TypedDataFrame[GumboWgsSequencing],
    task_results: TypedDataFrame[GumboTaskResult],
    ref_urls: dict[str, dict[str, str]],
    join_task_results: bool = False,
) -> TypedDataFrame[TerraSample]:
    """
    Make a data frame to upload to Terra as the `sample` data table.

    :param wgs_sequencings: data frame of Gumbo `omics_sequencing` records
    :param task_results: data frame of Gumbo `task_result` records
    :param ref_urls: a nested dictionary of genomes and their reference file URLs
    :param join_task_results: whether to join canonical task results to upserted samples
    :return: a data frame for the `sample` data table in Terra
    """

    # set up genome ref URLs as data frames
    hg38_urls = pd.DataFrame([ref_urls["hg38"]])
    hg38_urls["ref"] = "hg38"

    hg19_urls = pd.DataFrame([ref_urls["hg19"]])
    hg19_urls["ref"] = "hg19"

    # start constructing samples data frame
    samples = wgs_sequencings.copy()

    samples = samples.rename(
        columns={
            "bai_filepath": "analysis_ready_bai",
            "bam_filepath": "analysis_ready_bam",
            "sequencing_id": "sample_id",
        }
    )

    # indicate whether GP-delivered alignmnent files are hg19 or hg38
    samples["delivery_ref"] = pd.NA
    samples.loc[samples["hg38_cram_filepath"].notna(), "delivery_ref"] = "hg38"
    samples.loc[
        samples["delivery_ref"].isna() & samples["hg19_bam_filepath"].notna(),
        "delivery_ref",
    ] = "hg19"

    is_hg38 = samples["delivery_ref"].eq("hg38")
    is_hg19 = samples["delivery_ref"].eq("hg19")

    dfs = []

    # set delivery_* columns for hg38 samples
    if is_hg38.any():
        samples_sub = (
            samples.loc[is_hg38]
            .drop(columns=["hg19_bam_filepath", "hg19_bai_filepath"])
            .rename(
                columns={
                    "hg38_cram_filepath": "delivery_cram_bam",
                    "hg38_crai_filepath": "delivery_crai_bai",
                }
            )
        )

        hg38_delivery_urls = hg38_urls.copy()
        hg38_delivery_urls.columns = "delivery_" + hg38_delivery_urls.columns

        samples_sub = samples_sub.merge(
            hg38_delivery_urls, on="delivery_ref", how="left"
        )
        dfs.append(samples_sub)

    # set delivery_* columns for hg19 samples
    if is_hg19.any():
        samples_sub = (
            samples.loc[is_hg19]
            .drop(columns=["hg38_cram_filepath", "hg38_crai_filepath"])
            .rename(
                columns={
                    "hg19_bam_filepath": "delivery_cram_bam",
                    "hg19_bai_filepath": "delivery_crai_bai",
                }
            )
        )

        hg19_delivery_urls = hg19_urls.copy()
        hg19_delivery_urls.columns = "delivery_" + hg19_delivery_urls.columns

        samples_sub = samples_sub.merge(
            hg19_delivery_urls, on="delivery_ref", how="left"
        )
        dfs.append(samples_sub)

    # recombine hg38/19 samples
    samples = pd.concat(dfs, ignore_index=True)

    samples["delivery_file_format"] = samples["delivery_cram_bam"].apply(
        lambda x: pathlib.Path(x).suffix[1:].upper()
    )

    # for now, hardcoding the desired reference (for preprocessing/alignment workflow)
    samples["ref"] = "hg38"
    samples = samples.merge(hg38_urls, on="ref", how="left")

    if join_task_results:
        raise NotImplementedError

    return type_data_frame(samples, TerraSample)


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


def delta_preprocess_wgs_samples(
    terra_workspace: TerraWorkspace, terra_workflow: TerraWorkflow
) -> None:
    """
    Check the Terra `sample` data table for samples without `preprocess_wgs_sample`
    output files and submit a run of that workflow on a sample set comprised of those
    new files.

    :param terra_workspace: a `TerraWorkspace` instance
    :param terra_workflow: a `TerraWorkflow` instance
    """

    samples = terra_workspace.get_entities("sample", TerraSample)
    samples = samples.loc[
        samples["analysis_ready_bam"].isna() | samples["analysis_ready_bai"].isna()
    ]

    if len(samples) == 0:
        echo("No new samples to preprocess")
        return

    sample_set_id = terra_workspace.create_sample_set(
        samples["sample_id"], suffix="preprocess_wgs_sample"
    )

    terra_workspace.submit_workflow_run(
        terra_workflow,
        entity=sample_set_id,
        etype="sample_set",
        expression="this.samples",
        use_callcache=True,
        use_reference_disks=False,
        memory_retry_multiplier=1.2,
    )


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

        o = x.copy()

        # replace the entity name with the more reliable `sample_id` workflow input
        o.terra_entity_name = o.terra_workflow_inputs["sample_id"]

        # do trivial conversion from `TaskResult` to `task_result_insert_input` type
        outputs.append(task_result_insert_input(**o.model_dump()))

    echo("Getting existing task entity records for sequencings")
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
        echo(f"Creating {len(missing_seq_ids)} missing task entity records")
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

    echo(f"Upserting {len(outputs)} task results into Gumbo")
    res = gumbo_client.insert_terra_sync(
        username=gumbo_client.username,
        terra_workspace_namespace=str(outputs[0].terra_workspace_namespace),
        terra_workspace_name=str(outputs[0].terra_workspace_name),
        task_results=outputs,
    )

    sync_id = res.insert_terra_sync.returning[0].id  # pyright: ignore
    echo(f"Created Terra sync record {sync_id}")
