import datetime
import os
import pathlib
import tempfile
from typing import Any, Callable, Iterable

import pandas as pd
import requests
from click import echo
from firecloud import api as firecloud_api

from gumbo_gql_client import GumboClient, task_result_bool_exp
from omics_wgs_pipeline.types import TypedDataFrame
from omics_wgs_pipeline.utils import expand_dict_columns, model_to_df, type_data_frame
from omics_wgs_pipeline.validators import (
    GumboTaskResult,
    GumboWgsSequencing,
    SampleToPreprocess,
)


def call_firecloud_api(func: Callable, *args: object, **kwargs: object) -> Any:
    """
    Call a Firecloud API endpoint and check the response for a valid HTTP status code.

    :param func: a `firecloud.api` method
    :param args: arguments to `func`
    :param kwargs: keyword arguments to `func`
    :return: the API response, if any
    """

    res = func(*args, **kwargs)

    if 200 <= res.status_code <= 299:
        try:
            return res.json()
        except requests.exceptions.JSONDecodeError:
            return res.text

    try:
        raise requests.exceptions.RequestException(
            f"HTTP {res.status_code} error: {res.json()}"
        )
    except Exception as e:
        # it's returning HTML or something is preventing us from parsing the JSON
        echo(f"Error getting response as JSON: {e}", err=True)
        echo(f"Response text: {res.text}", err=True)
        raise requests.exceptions.RequestException(f"HTTP {res.status_code} error")


def upload_entities_to_terra(
    workspace_namespace: str, workspace_name: str, df: pd.DataFrame
) -> None:
    """
    Upload a data frame of entities to a Terra data table.

    :param workspace_namespace: the namespace of the Firecloud workspace
    :param workspace_name: the name of the Firecloud workspace
    :param df: a data frame of entities
    """

    with tempfile.NamedTemporaryFile(suffix="tsv") as f:
        df.to_csv(f, sep="\t", index=False)

        call_firecloud_api(
            firecloud_api.upload_entities_tsv,
            namespace=workspace_namespace,
            workspace=workspace_name,
            entities_tsv=f.name,
            model="flexible",
        )


def refresh_terra_samples(workspace_namespace: str, workspace_name: str) -> None:
    """
    Update the sample data table using ground truth data from Gumbo.

    :param workspace_namespace: the namespace of the Firecloud workspace
    :param workspace_name: the name of the Firecloud workspace
    """

    # collect ground truth set of samples from Gumbo
    gumbo_client = GumboClient(
        url=os.environ["HASURA_URL"],
        headers={"X-Hasura-Admin-Secret": os.environ["HASURA_ADMIN_SECRET"]},
    )

    wgs_sequencings = model_to_df(
        gumbo_client.wgs_sequencings(),
        GumboWgsSequencing,
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
            task_result_bool_exp(workflow_name={"eq": "preprocess_wgs_sample"})
        ),
        GumboTaskResult,
        mutator=lambda df: expand_dict_columns(df).rename(
            columns={"task_entity__sequencing_id": "sample_id"}
        ),
    )

    samples = make_samples_to_preprocess(wgs_sequencings, task_results)

    echo(f"Upserting {len(samples)} samples to Terra")
    upload_entities_to_terra(workspace_namespace, workspace_name, df=samples)


def make_samples_to_preprocess(
    wgs_sequencings: TypedDataFrame[GumboWgsSequencing],
    task_results: TypedDataFrame[GumboTaskResult],
) -> TypedDataFrame[SampleToPreprocess]:
    """
    Make a data frame to upload to Terra as the `sample` data table.

    :param wgs_sequencings: data frame of Gumbo `omics_sequencing` records
    :param task_results: data frame of Gumbo `task_result` records
    :return: a data frame for the `sample` data table in Terra
    """

    samples = wgs_sequencings.copy()

    samples = samples.rename(columns={"sequencing_id": "sample_id"})

    # pick the best available option for the alignment files
    samples["delivery_cram_bam"] = (
        samples["hg38_cram_filepath"]
        .fillna(samples["bam_filepath"])
        .fillna(samples["hg19_bam_filepath"])
    )

    samples["delivery_crai_bai"] = (
        samples["hg38_crai_filepath"]
        .fillna(samples["bai_filepath"])
        .fillna(samples["hg19_bai_filepath"])
    )

    samples["delivery_file_format"] = samples["delivery_cram_bam"].apply(
        lambda x: pathlib.Path(x).suffix[1:].upper()
    )

    # all new samples should be hg38, but this removes hg19 ones if necessary
    samples["is_hg19"] = (
        samples["delivery_cram_bam"].eq(samples["hg19_bam_filepath"]).fillna(False)
    )
    samples = samples.loc[~samples["is_hg19"]]

    # join already preprocessed BAM/BAI files to canonical set of WGS samples
    bam_bais = task_results.pivot(
        index="sample_id", columns="label", values="url"
    ).reset_index()

    samples = samples.merge(bam_bais, how="left", on="sample_id")

    return type_data_frame(samples, SampleToPreprocess)


def create_new_sample_set(
    workspace_namespace: str,
    workspace_name: str,
    sample_ids: Iterable[str],
    suffix: str,
) -> str:
    """
    Create a new sample set for samples and upload it to Terra.

    :param workspace_namespace: the namespace of the Firecloud workspace
    :param workspace_name: the name of the Firecloud workspace
    :param sample_ids: a list of sample IDs
    :param suffix: a suffix to add to the sample set ID (e.g. "preprocess_wgs_sample")
    :return: the ID of the new sample set
    """

    # make an ID for the sample set of new samples
    sample_set_id = "_".join(
        [
            "samples",
            datetime.datetime.now(datetime.UTC)
            .isoformat(timespec="seconds")
            .replace(":", "-"),
            suffix,
        ]
    )

    # construct a data frame of sample IDs for this sample set
    sample_sets = pd.DataFrame({"entity:sample_id": sample_ids}, dtype="string")
    sample_sets["entity:sample_set_id"] = sample_set_id

    echo("Creating new sample set in Terra")
    upload_entities_to_terra(
        workspace_namespace,
        workspace_name,
        df=sample_sets[["entity:sample_set_id"]].drop_duplicates(),
    )

    # construct the join table between the sample set and its samples
    sample_sets = sample_sets.rename(
        columns={
            "entity:sample_set_id": "membership:sample_set_id",
            "entity:sample_id": "sample",
        }
    )

    sample_sets = sample_sets.loc[:, ["membership:sample_set_id", "sample"]]

    echo(f"Adding {len(sample_sets)} samples to sample set {sample_set_id}")
    upload_entities_to_terra(
        workspace_namespace,
        workspace_name,
        df=sample_sets,
    )

    return sample_set_id


def submit_workflow_run(
    workspace_namespace: str,
    workspace_name: str,
    repo_namespace: str,
    repo_method_name: str,
    sample_set_id: str,
) -> None:
    """
    Submit a run of a workflow for the new sample set.

    :param workspace_namespace: the namespace of the Firecloud workspace
    :param workspace_name: the name of the Firecloud workspace
    :param repo_namespace: the namespace of the Firecloud method
    :param repo_method_name: the name of the Firecloud method
    :param sample_set_id: the ID of the new sample set
    """

    call_firecloud_api(
        firecloud_api.create_submission,
        wnamespace=workspace_namespace,
        workspace=workspace_name,
        cnamespace=repo_namespace,
        config=repo_method_name,
        entity=sample_set_id,
        etype="sample_set",
        expression="this.samples",
        use_callcache=True,
        use_reference_disks=True,
        memory_retry_multiplier=1.2,  # pyright: ignore
    )
