import logging
from functools import partial

import pandas as pd
import requests
from firecloud_api_cds import api as firecloud_api
from nebelung.terra_workspace import TerraWorkspace
from nebelung.utils import call_firecloud_api, type_data_frame
from pd_flatten import pd_flatten

from depmap_omics_wgs.gcs import copy_to_cclebams, delete_blobs, get_objects_metadata
from depmap_omics_wgs.types import (
    AlignedSamplesWithObjectMetadata,
    GumboClient,
    GumboWgsSequencing,
    NewSequencingAlignment,
    TerraSample,
    TypedDataFrame,
    UpdatedSequencingAlignment,
)
from depmap_omics_wgs.utils import model_to_df
from gumbo_gql_client import (
    sequencing_alignment_insert_input,
    sequencing_alignment_set_input,
)


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
        .drop(columns=["sequencing_alignment_source"])
        .rename(
            columns={
                "omics_sequencing_id": "sample_id",
                "sequencing_alignment_id": "delivery_sequencing_alignment_id",
                "url": "delivery_cram_bam",
                "size": "delivery_cram_bam_size",
                "index_url": "delivery_crai_bai",
                "reference_genome": "delivery_ref",
            }
        )
        .merge(
            pd.concat(
                [
                    wgs_sequencings.loc[
                        wgs_sequencings["sequencing_alignment_source"].eq("CDS")
                        & wgs_sequencings["url"].str.endswith(".bam")
                    ]
                    .drop(
                        columns=[
                            "sequencing_alignment_source",
                            "model_id",
                            "model_condition_id",
                            "omics_profile_id",
                            "cell_line_name",
                            "size",
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
                    wgs_sequencings.loc[
                        wgs_sequencings["sequencing_alignment_source"].eq("CDS")
                        & wgs_sequencings["url"].str.endswith(".cram")
                    ]
                    .drop(
                        columns=[
                            "sequencing_alignment_source",
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
                            "url": "analysis_ready_cram",
                            "size": "analysis_ready_cram_size",
                            "index_url": "analysis_ready_crai",
                            "reference_genome": "ref",
                        }
                    ),
                ],
                ignore_index=True,
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

    # only newly loaded samples can be processed by automation
    samples["automation_status"] = pd.NA

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
    # don't replace recently populated cells like `analysis_ready_cram` with blanks in
    # case we haven't run `onboard_aligned_crams` since the most recent alignment jobs
    # have completed, thus `delete_empty=False`
    terra_workspace.upload_entities(df=samples, delete_empty=False)

    # if a sample has an analysis-ready CRAM, though, that means we've completed all
    # workflows, converted the temporary analysis-ready BAM to CRAM, and archived the
    # CRAM, so we can blank out the analysis_ready_ba{i,m} columns now
    samples_done = samples.dropna(subset=["analysis_ready_cram"]).loc[
        :, ["entity:sample_id", "analysis_ready_bai", "analysis_ready_bam"]
    ]
    terra_workspace.upload_entities(df=samples_done, delete_empty=True)


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
) -> None:
    """
    Update the sample data table in the legacy production Terra workspace.

    :param terra_workspace: the Terra workspace where new workflows have run
    :param legacy_terra_workspace: the legacy Terra workspace where new workflow outputs
    need to be synced to
    :param sample_set_id: the ID for a sample set to create or update for unprocessed
    samples
    """

    sample_sets = terra_workspace.get_entities("sample_set")
    sample_set = sample_sets.loc[
        sample_sets["sample_set_id"].eq(sample_set_id), "samples"
    ].values[0]
    q_sample_ids = [x["entityName"] for x in sample_set["items"]]

    src_samples = terra_workspace.get_entities("sample")

    dest_samples = (
        src_samples.loc[
            :,
            [
                "sample_id",
                "cnv_cn_by_gene_weighted_mean",
                "cnv_segments",
                "mut_somatic_variants",
                "sv_selected_somatic",
                "msisensor2_output_dis",
                "msisensor2_score",
                "guide_bed_avana",
                "guide_bed_brunello",
                "guide_bed_humagne",
                "guide_bed_ky",
                "guide_bed_tkov",
                "mut_annot_bcftools_fixed_vcf",
            ],
        ]
        .rename(
            columns={
                "mut_somatic_variants": "depmap_maf_25q2",
                "sv_selected_somatic": "expanded_filtered_sv_bedpe",
                "guide_bed_avana": "avana_binary_mut",
                "guide_bed_brunello": "brunello_binary_mut",
                "guide_bed_humagne": "humagne_binary_mut",
                "guide_bed_ky": "ky_binary_mut",
                "guide_bed_tkov": "tkov3_binary_mut",
                "mut_annot_bcftools_fixed_vcf": "mutect2_fixed_vcf",
            }
        )
        .astype("string")
        .replace({"": pd.NA})
    )

    dest_samples = dest_samples.loc[
        dest_samples.drop(columns="sample_id").notna().any(axis=1)
    ]

    if len(dest_samples) > 0:
        legacy_terra_workspace.upload_entities(dest_samples, delete_empty=False)

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
            {"entityType": "sample", "entityName": x} for x in q_sample_ids
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
            entity_type="sample", entity_ids=q_sample_ids, entity_set_id=sample_set_id
        )

    try:
        # update existing sample set
        sample_set = call_firecloud_api(
            firecloud_api.get_entity,
            namespace=legacy_terra_workspace.workspace_namespace,
            workspace=legacy_terra_workspace.workspace_name,
            etype="sample_set",
            ename="all_25q3",
        )

        sample_set["attributes"]["samples"]["items"] = [
            {"entityType": "sample", "entityName": x} for x in dest_samples["sample_id"]
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
            entity_ids=list(dest_samples["sample_id"]),
            entity_set_id="all_25q3",
        )


def onboard_aligned_crams(
    gumbo_client: GumboClient,
    terra_workspace: TerraWorkspace,
    gcp_project_id: str,
    ref_urls: dict[str, dict[str, str]],
    dry_run: bool,
) -> None:
    """
    Copy analysis-ready CRAM/CRAI files to cclebams and onboard sequencing alignment
    records to Gumbo.

    :param terra_workspace: a TerraWorkspace instance
    :param gumbo_client: an instance of the Gumbo GraphQL client
    :param gcp_project_id: the ID of a GCP project to use for billing
    :param ref_urls: a nested dictionary of genomes and their reference file URLs
    :param dry_run: whether to skip updates to external data stores
    """

    # get the alignment files from Terra
    samples = terra_workspace.get_entities("sample")

    # subset to samples with all workflows completed (i.e. we're ready to archive the
    # analysis CRAM now)
    samples = (
        samples.dropna(
            subset=[
                "sample_id",
                "analysis_ready_crai",
                "analysis_ready_cram",
                "cnv_bin_coverage",
                "cnv_cn_by_gene_weighted_mean",
                "cnv_input_params",
                "cnv_read_cov_bin",
                "cnv_segments",
                "guide_bed_avana",
                "guide_bed_brunello",
                "guide_bed_humagne",
                "guide_bed_ky",
                "guide_bed_tkov",
                "msisensor2_output",
                "msisensor2_output_dis",
                "msisensor2_output_somatic",
                "msisensor2_score",
                "mut_annot_bcftools_fixed_vcf",
                "mut_annot_bcftools_vcf",
                "mut_annot_open_cravat_vcf",
                "mut_annot_snpeff_snpsift_vcf",
                "mut_annot_vcf",
                "mut_annot_vcf_index",
                "mut_annot_vep_vcf",
                "mut_duckdb",
                "mut_filtering_stats",
                "mut_sig_variants",
                "mut_somatic_variants",
                "mut_vcf",
                "mut_vcf_idx",
                "sv_annot_bedpe",
                "sv_annot_reannotated_bedpe",
                "sv_annot_vcf",
                "sv_candidate_indel_vcf",
                "sv_candidate_indel_vcf_index",
                "sv_candidate_vcf",
                "sv_candidate_vcf_index",
                "sv_del_annotation",
                "sv_dup_annotation",
                "sv_selected_somatic",
                "sv_somatic_vcf",
                "sv_somatic_vcf_index",
            ]
        )
        .loc[
            :,
            [
                "sample_id",
                "analysis_ready_bam",
                "analysis_ready_bai",
                "analysis_ready_cram",
                "analysis_ready_crai",
            ],
        ]
        .rename(columns={"sample_id": "omics_sequencing_id"})
    )

    # get sequencing alignment records for WGS omics_sequencings from Gumbo
    existing_alignments = model_to_df(
        gumbo_client.wgs_sequencing_alignments(timeout=30.0),
        GumboWgsSequencing,
        mutator=partial(pd_flatten, name_columns_with_parent=False),
    )

    # subset to the analysis-ready records
    existing_alignments = existing_alignments.loc[
        existing_alignments["sequencing_alignment_source"].eq("CDS"),
        ["omics_sequencing_id", "sequencing_alignment_id", "url", "index_url", "size"],
    ]

    # compare Gumbo vs. Terra
    comp = existing_alignments.merge(
        samples, how="outer", on="omics_sequencing_id"
    ).drop(columns="size")

    # identify records in Gumbo that are missing or need to be updated (e.g. replacing
    # an analysis-ready BAM with a CRAM)
    comp = comp.loc[comp["url"].isna() | comp["url"].ne(comp["analysis_ready_cram"])]

    if len(comp) == 0:
        logging.info("No new/changed aligned CRAM files to onboard")
        return

    # get GCS blob metadata for the CRAMs
    objects_metadata = get_objects_metadata(comp["analysis_ready_cram"])

    comp = type_data_frame(
        comp.merge(
            objects_metadata.rename(columns={"url": "analysis_ready_cram"}),
            how="left",
            on="analysis_ready_cram",
        ).drop(columns=["gcs_obj_updated_at"]),
        AlignedSamplesWithObjectMetadata,
    )

    # confirm again using file size that the CRAMs don't already exist as Gumbo records
    assert not bool(comp["size"].isin(existing_alignments["size"]).any())

    # copy CRAMs and CRAIs to our bucket
    copied_files = copy_to_cclebams(
        comp,
        gcp_project_id=gcp_project_id,
        gcs_destination_bucket="cclebams",
        gcs_destination_prefix="wgs_hg38",
        dry_run=dry_run,
    )

    assert bool(copied_files["copied"].all())

    # delete obsolete alignment files
    blobs_to_delete = (
        pd.concat(
            [
                # can delete the workspace's analysis ready BAM
                comp["analysis_ready_bam"],
                comp["analysis_ready_bai"],
                # can delete the old analysis ready BAM/CRAM referenced in Gumbo
                comp["url"],
                comp["index_url"],
            ]
        )
        .dropna()
        .drop_duplicates()
    )

    if len(blobs_to_delete) > 0:
        delete_blobs(blobs_to_delete, gcp_project_id, dry_run)

    # update URLs with the ones we just copied to the archive bucket
    samples_updated = comp.copy().drop(columns=["url", "index_url"])

    for c in ["analysis_ready_crai", "analysis_ready_cram"]:
        sample_file_urls = copied_files.loc[copied_files["copied"], ["url", "new_url"]]
        samples_updated = samples_updated.merge(
            sample_file_urls, how="left", left_on=c, right_on="url"
        )
        samples_updated[c] = samples_updated["new_url"]
        samples_updated = samples_updated.drop(columns=["url", "new_url"])

    # create/update sequencing_alignment records in Gumbo
    persist_sequencing_alignments(
        gumbo_client, samples=samples_updated, dry_run=dry_run
    )

    # ensure all BAM/CRAM cols in workspace are up to date
    refresh_terra_samples(
        terra_workspace=terra_workspace,
        gumbo_client=gumbo_client,
        ref_urls=config["ref"],
    )


def persist_sequencing_alignments(
    gumbo_client: GumboClient,
    samples: TypedDataFrame[AlignedSamplesWithObjectMetadata],
    dry_run: bool,
) -> None:
    """
    Insert many sequencing alignments to Gumbo.

    :param gumbo_client: an instance of the Gumbo GraphQL client
    :param samples: a data frame of prepared sequencing alignments to insert/update
    :param dry_run: whether to skip updates to external data stores
    """

    # set columns and names for Gumbo inserts/updates
    samples["reference_genome"] = "hg38"
    samples["sequencing_alignment_source"] = "CDS"

    samples = samples.rename(
        columns={
            "crc32c": "crc32c_hash",
            "analysis_ready_cram": "url",
            "analysis_ready_crai": "index_url",
        }
    )

    # collect and insert new records
    new_sequencing_alignments = type_data_frame(
        samples.loc[samples["sequencing_alignment_id"].isna()],
        NewSequencingAlignment,
        remove_unknown_cols=True,
    )

    if len(new_sequencing_alignments) > 0:
        persist_new_sequencing_alignments(
            gumbo_client, samples=new_sequencing_alignments, dry_run=dry_run
        )

    # collect and update changed records
    updated_sequencing_alignments = type_data_frame(
        samples.loc[samples["sequencing_alignment_id"].notna()],
        UpdatedSequencingAlignment,
    )

    if len(updated_sequencing_alignments) > 0:
        persist_updated_sequencing_alignments(
            gumbo_client, samples=updated_sequencing_alignments, dry_run=dry_run
        )


def persist_new_sequencing_alignments(
    gumbo_client: GumboClient,
    samples: TypedDataFrame[NewSequencingAlignment],
    dry_run: bool = True,
):
    """
    Insert many sequencing alignments to Gumbo.

    :param gumbo_client: an instance of the Gumbo GraphQL client
    :param samples: a data frame of prepared sequencing alignments to insert
    :param dry_run: whether to skip updates to external data stores
    """

    if dry_run:
        logging.info(f"(skipping) Inserting {len(samples)} sequencing alignments")
        return

    sequencing_alignment_inserts = [
        sequencing_alignment_insert_input.model_validate(x)
        for x in samples.to_dict(orient="records")
    ]

    res = gumbo_client.insert_sequencing_alignments(
        gumbo_client.username, objects=sequencing_alignment_inserts, timeout=60.0
    )

    affected_rows = res.insert_sequencing_alignment.affected_rows  # pyright: ignore
    logging.info(f"Inserted {affected_rows} sequencing alignments")


def persist_updated_sequencing_alignments(
    gumbo_client: GumboClient,
    samples: TypedDataFrame[UpdatedSequencingAlignment],
    dry_run: bool = True,
):
    """
    Updated many sequencing alignments in Gumbo.

    :param gumbo_client: an instance of the Gumbo GraphQL client
    :param samples: a data frame of prepared sequencing alignments to update
    :param dry_run: whether to skip updates to external data stores
    """

    if dry_run:
        logging.info(f"(skipping) Updating {len(samples)} sequencing alignments")
        return

    for _, r in samples.iterrows():
        sequencing_alignment_update = sequencing_alignment_set_input.model_validate(
            r.to_dict()
        )

        # noinspection PyTypeChecker
        res = gumbo_client.update_sequencing_alignment(
            gumbo_client.username,
            id=r["sequencing_alignment_id"],
            object=sequencing_alignment_update,
            timeout=30.0,
        )

        updated_id = res.update_sequencing_alignment_by_pk.id  # pyright: ignore
        logging.info(f"Updated sequencing alignment {updated_id}")
