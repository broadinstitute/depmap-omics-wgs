import logging
from functools import partial

import pandas as pd
import requests
from firecloud_api_cds import api as firecloud_api
from nebelung.terra_workspace import TerraWorkspace
from nebelung.utils import call_firecloud_api, type_data_frame
from pd_flatten import pd_flatten

from depmap_omics_wgs.gcs import copy_to_cclebams, get_objects_metadata
from depmap_omics_wgs.types import (
    AlignedSamplesWithObjectMetadata,
    CopiedSampleFiles,
    GumboClient,
    GumboWgsSequencing,
    NewSequencingAlignments,
    TerraSample,
    TypedDataFrame,
)
from depmap_omics_wgs.utils import model_to_df
from gumbo_gql_client import sequencing_alignment_insert_input


def refresh_terra_samples(
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
        gumbo_client.wgs_sequencing_alignments(timeout=30.0),
        GumboWgsSequencing,
        mutator=partial(pd_flatten, name_columns_with_parent=False),
    )

    # make wide, separating delivery and analysis-ready CRAM/BAMs
    samples = (
        wgs_sequencings.loc[wgs_sequencings["sequencing_alignment_source"].eq("GP")]
        .drop(columns=["sequencing_alignment_source", "size"])
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
            .drop(
                columns=[
                    "sequencing_alignment_source",
                    "size",
                    "model_id",
                    "model_condition_id",
                    "omics_profile_id",
                    "cell_line_name",
                    "stripped_cell_line_name",
                ]
            )
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
    # don't replace recently populated cells like `analysis_ready_bam` with blanks in
    # case we haven't run `onboard_aligned_bams` since the most recent alignment jobs
    # have completed, thus `delete_empty=False`
    terra_workspace.upload_entities(df=samples, delete_empty=False)


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


def refresh_legacy_terra_samples(
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

    src_samples = (
        samples.loc[
            :,
            [
                "sample_id",
                "analysis_ready_bam",
                "analysis_ready_bai",
                "cnv_cn_by_gene_weighted_mean",
                "cnv_segments",
            ],
        ]
        .rename(
            columns={
                "analysis_ready_bam": "internal_bam_filepath",
                "analysis_ready_bai": "internal_bai_filepath",
            }
        )
        .replace({"": pd.NA})
    )

    if "internal_bam_filepath" not in legacy_samples.columns:
        legacy_samples["internal_bam_filepath"] = pd.NA

    if "internal_bai_filepath" not in legacy_samples.columns:
        legacy_samples["internal_bai_filepath"] = pd.NA

    if "cnv_cn_by_gene_weighted_mean" not in legacy_samples.columns:
        legacy_samples["cnv_cn_by_gene_weighted_mean"] = pd.NA

    if "cnv_segments" not in legacy_samples.columns:
        legacy_samples["cnv_segments"] = pd.NA

    dest_samples = legacy_samples.loc[
        :,
        [
            "sample_id",
            "internal_bam_filepath",
            "internal_bai_filepath",
            "cnv_cn_by_gene_weighted_mean",
            "cnv_segments",
        ],
    ].replace({"": pd.NA})

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
    res = gumbo_client.unprocessed_sequencings(timeout=30.0)

    unprocessed_sample_ids = [
        x.id for x in res.records if x.id in legacy_samples["sample_id"].values
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


def onboard_aligned_bams(
    gumbo_client: GumboClient,
    terra_workspace: TerraWorkspace,
    gcp_project_id: str,
    dry_run: bool,
) -> None:
    """
    Copy BAM/BAI files to cclebams and onboard sequencing alignment records to Gumbo.

    :param terra_workspace: a TerraWorkspace instance
    :param gumbo_client: an instance of the Gumbo GraphQL client
    :param gcp_project_id: the ID of a GCP project to use for billing
    :param dry_run: whether to skip updates to external data stores
    """

    # get the alignment files from Terra
    samples = terra_workspace.get_entities("sample")
    samples = samples.loc[
        :, ["sample_id", "analysis_ready_bam", "analysis_ready_bai"]
    ].dropna()

    # get sequencing alignment records for WGS omics_sequencings from Gumbo
    existing_alignments = model_to_df(
        gumbo_client.wgs_sequencing_alignments(timeout=30.0),
        GumboWgsSequencing,
        mutator=partial(pd_flatten, name_columns_with_parent=False),
    )

    # check which sequencings have delivery BAMs but not analysis-ready/aligned BAMs
    seq_ids_no_cds = set(
        existing_alignments.loc[
            existing_alignments["sequencing_alignment_source"].eq("GP"),
            "omics_sequencing_id",
        ]
    ).difference(
        existing_alignments.loc[
            existing_alignments["sequencing_alignment_source"].eq("CDS"),
            "omics_sequencing_id",
        ]
    )

    if len(seq_ids_no_cds) == 0:
        logging.info("No new aligned BAM files to onboard")
        return

    # subset to newly-aligned sequencings that need to be onboarded
    samples = samples.loc[samples["sample_id"].isin(list(seq_ids_no_cds))]

    # get GCS blob metadata for the BAMs
    objects_metadata = get_objects_metadata(samples["analysis_ready_bam"])

    samples = type_data_frame(
        samples.merge(
            objects_metadata, how="left", left_on="analysis_ready_bam", right_on="url"
        ).drop(columns=["url", "gcs_obj_updated_at"]),
        AlignedSamplesWithObjectMetadata,
    )

    # confirm again using file size that these BAMs don't already exist as Gumbo records
    assert ~bool(samples["size"].isin(existing_alignments["size"]).any())

    # copy BAMs and BAIs to our bucket
    sample_files = copy_to_cclebams(
        samples,
        gcp_project_id=gcp_project_id,
        gcs_destination_bucket="cclebams",
        gcs_destination_prefix="wgs_hg38",
        dry_run=dry_run,
    )

    samples = update_sample_file_urls(samples, sample_files)

    # create sequencing_alignment records in Gumbo
    persist_sequencing_alignments(gumbo_client, samples, dry_run)


def update_sample_file_urls(
    samples: TypedDataFrame[AlignedSamplesWithObjectMetadata],
    sample_files: TypedDataFrame[CopiedSampleFiles],
) -> TypedDataFrame[AlignedSamplesWithObjectMetadata]:
    """
    Replace BAM URLs with new ones used in `copy_to_depmap_omics_bucket`.

    :param samples: the data frame of samples
    :param sample_files: a data frame of files we attempted to copy
    :return: the samples data frame with issue column filled out for rows with files we
    couldn't copy
    """

    logging.info("Updating GCS file URLs...")

    samples_updated = samples.copy()

    for c in ["analysis_ready_bai", "analysis_ready_bam"]:
        sample_file_urls = sample_files.loc[sample_files["copied"], ["url", "new_url"]]
        samples_updated = samples_updated.merge(
            sample_file_urls, how="left", left_on=c, right_on="url"
        )
        samples_updated[c] = samples_updated["new_url"]
        samples_updated = samples_updated.drop(columns=["url", "new_url"])

    return type_data_frame(samples_updated, AlignedSamplesWithObjectMetadata)


def persist_sequencing_alignments(
    gumbo_client: GumboClient,
    samples: TypedDataFrame[AlignedSamplesWithObjectMetadata],
    dry_run: bool,
) -> None:
    """
    Insert many sequencing alignments to Gumbo.

    :param gumbo_client: an instance of the Gumbo GraphQL client
    :param samples: a data frame of prepared sequencing alignments to insert
    :param dry_run: whether to skip updates to external data stores
    """

    new_sequencing_alignments = samples.rename(
        columns={
            "sample_id": "omics_sequencing_id",
            "crc32c": "crc32c_hash",
            "analysis_ready_bam": "url",
            "analysis_ready_bai": "index_url",
        }
    )

    new_sequencing_alignments["reference_genome"] = "hg38"
    new_sequencing_alignments["sequencing_alignment_source"] = "CDS"

    new_sequencing_alignments = type_data_frame(
        new_sequencing_alignments, NewSequencingAlignments
    )

    sequencing_alignment_inserts = [
        sequencing_alignment_insert_input.model_validate(x)
        for x in new_sequencing_alignments.to_dict(orient="records")
    ]

    if dry_run:
        logging.info(
            f"(skipping) Inserting {len(sequencing_alignment_inserts)} "
            "sequencing alignments"
        )
        return

    res = gumbo_client.insert_sequencing_alignments(
        gumbo_client.username, objects=sequencing_alignment_inserts, timeout=60.0
    )

    affected_rows = res.insert_sequencing_alignment.affected_rows  # pyright: ignore
    logging.info(f"Inserted {affected_rows} sequencing alignments")
