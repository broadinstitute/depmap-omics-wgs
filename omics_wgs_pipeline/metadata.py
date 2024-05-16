import os
import pathlib

from gumbo_gql_client import GumboClient, task_result_bool_exp
from omics_wgs_pipeline.types import TypedDataFrame
from omics_wgs_pipeline.utils import expand_dict_columns, model_to_df, type_data_frame
from omics_wgs_pipeline.validators import (
    GumboTaskResult,
    GumboWgsSequencing,
    TerraSample,
)


def make_terra_samples() -> TypedDataFrame[TerraSample]:
    """
    Make a data frame to use as a Terra `sample` data table using ground truth data from
    Gumbo.

    :return: a data frame to use as a Terra `sample` data table
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
            # TODO: other output files
            task_result_bool_exp(workflow_name={"eq": "preprocess_wgs_sample"})
        ),
        GumboTaskResult,
        mutator=lambda df: expand_dict_columns(df).rename(
            columns={"task_entity__sequencing_id": "sample_id"}
        ),
    )

    return join_existing_results_to_samples(wgs_sequencings, task_results)


def join_existing_results_to_samples(
    wgs_sequencings: TypedDataFrame[GumboWgsSequencing],
    task_results: TypedDataFrame[GumboTaskResult],
) -> TypedDataFrame[TerraSample]:
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

    # join already processed output files to canonical set of WGS samples
    task_result_urls = task_results.pivot(
        index="sample_id", columns="label", values="url"
    ).reset_index()

    samples = samples.merge(task_result_urls, how="left", on="sample_id")

    return type_data_frame(samples, TerraSample)


# def create_new_sample_set(
#     workspace_namespace: str,
#     workspace_name: str,
#     sample_ids: Iterable[str],
#     suffix: str,
# ) -> str:
#     """
#     Create a new sample set for samples and upload it to Terra.
#
#     :param workspace_namespace: the namespace of the Firecloud workspace
#     :param workspace_name: the name of the Firecloud workspace
#     :param sample_ids: a list of sample IDs
#     :param suffix: a suffix to add to the sample set ID (e.g. "preprocess_wgs_sample")
#     :return: the ID of the new sample set
#     """
#
#     # make an ID for the sample set of new samples
#     sample_set_id = "_".join(
#         [
#             "samples",
#             datetime.datetime.now(datetime.UTC)
#             .isoformat(timespec="seconds")
#             .replace(":", "-"),
#             suffix,
#         ]
#     )
#
#     # construct a data frame of sample IDs for this sample set
#     sample_sets = pd.DataFrame({"entity:sample_id": sample_ids}, dtype="string")
#     sample_sets["entity:sample_set_id"] = sample_set_id
#
#     echo("Creating new sample set in Terra")
#     upload_entities_to_terra(
#         workspace_namespace,
#         workspace_name,
#         df=sample_sets[["entity:sample_set_id"]].drop_duplicates(),
#     )
#
#     # construct the join table between the sample set and its samples
#     sample_sets = sample_sets.rename(
#         columns={
#             "entity:sample_set_id": "membership:sample_set_id",
#             "entity:sample_id": "sample",
#         }
#     )
#
#     sample_sets = sample_sets.loc[:, ["membership:sample_set_id", "sample"]]
#
#     echo(f"Adding {len(sample_sets)} samples to sample set {sample_set_id}")
#     upload_entities_to_terra(workspace_namespace, workspace_name, df=sample_sets)
#
#     return sample_set_id
