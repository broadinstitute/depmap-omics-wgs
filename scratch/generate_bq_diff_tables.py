from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd
from google.cloud import bigquery, storage
from nebelung.terra_workspace import TerraWorkspace

pd.set_option("display.max_columns", 30)
pd.set_option("display.max_colwidth", 50)
pd.set_option("display.max_info_columns", 30)
pd.set_option("display.max_info_rows", 20)
pd.set_option("display.max_rows", 20)
pd.set_option("display.max_seq_items", None)
pd.set_option("display.width", 200)
pd.set_option("expand_frame_repr", True)
pd.set_option("mode.chained_assignment", "warn")

ws = TerraWorkspace(
    workspace_namespace="broad-firecloud-ccle",
    workspace_name="depmap-omics-wgs-mut-dev",
)

samples = ws.get_entities("sample")

client = bigquery.Client()
dataset_id = "depmap-omics.mut_eval"


for sample_id in samples["sample_id"]:
    table_suffix = sample_id.replace("-", "_")
    old_table = f"{dataset_id}.{table_suffix}_somatic"

    query = f"""
        WITH
          new_muts AS (
            SELECT * FROM `{dataset_id}.new`
            WHERE sample_id = '{sample_id}' AND somatic
          ),
          old_muts AS (
            SELECT chrom, pos, ref, alt FROM `{old_table}`
          )
        SELECT
          *
        FROM
          new_muts
        LEFT JOIN
          old_muts
        USING (chrom, pos, ref, alt)
        WHERE
          old_muts.chrom IS NULL
    """

    destination_table = f"{dataset_id}.only_new_new"

    job_config = bigquery.QueryJobConfig(
        destination=destination_table,
        write_disposition="WRITE_APPEND",
    )

    print(f"Inserting into only_new_new for {sample_id}")
    query_job = client.query(query, job_config=job_config)
    query_job.result()

for sample_id in samples["sample_id"]:
    table_suffix = sample_id.replace("-", "_")
    old_table = f"{dataset_id}.{table_suffix}_somatic"

    query = f"""
        WITH
          new_muts AS (
          SELECT
            chrom,
            pos,
            ref,
            alt
          FROM
            `{dataset_id}.new`
          WHERE
            sample_id='{sample_id}'
            AND somatic ),
          old_muts AS (
          SELECT
            cast(chrom AS STRING) as chrom,
            cast(pos AS INT64) as pos,
            cast(ref AS STRING) as ref,
            cast(alt AS STRING) as alt,
            cast(af AS FLOAT64) as af,
            cast(dp AS INT64) as dp,
            cast(ref_count AS INT64) as ref_count,
            cast(alt_count AS INT64) as alt_count,
            cast(gt AS STRING) as gt,
            cast(ps AS INT64) as ps,
            cast(variant_type AS STRING) as variant_type,
            cast(variant_info AS STRING) as variant_info,
            cast(dna_change AS STRING) as dna_change,
            cast(protein_change AS STRING) as protein_change,
            cast(hugo_symbol AS STRING) as hugo_symbol,
            cast(hgnc_name AS STRING) as hgnc_name,
            cast(hgnc_family AS STRING) as hgnc_family,
            cast(vep_impact AS STRING) as vep_impact,
            cast(gnomade_af AS FLOAT64) as gnomade_af,
            cast(gnomadg_af AS FLOAT64) as gnomadg_af,
            cast(oncogene_high_impact AS BOOL) as oncogene_high_impact,
            cast(tumor_suppressor_high_impact AS BOOL) as tumor_suppressor_high_impact,
            cast(lof AS BOOL) as lof,
            cast(driver AS BOOL) as driver,
            cast(likely_driver AS BOOL) as likely_driver,
            cast(transcript_likely_lof AS STRING) as transcript_likely_lof,
            cast(brca1_func_score AS FLOAT64) as brca1_func_score,
            cast(civic_id AS FLOAT64) as civic_id,
            cast(civic_description AS STRING) as civic_description,
            cast(civic_score AS FLOAT64) as civic_score,
            cast(hess_driver AS BOOL) as hess_driver,
            cast(hess_signture AS STRING) as hess_signture,
            cast(cosmic_tier AS INT64) as cosmic_tier,
            cast(oncokb_effect AS STRING) as oncokb_effect,
            cast(oncokb_hotspot AS BOOL) as oncokb_hotspot,
            cast(oncokb_oncogenic AS STRING) as oncokb_oncogenic,
            cast(segdup AS BOOL) as segdup,
            cast(rm AS BOOL) as rm,
            cast(rescue AS BOOL) as rescue
          FROM
            `{old_table}` )
        SELECT
          '{sample_id}' as sample_id,
          *
        FROM
          old_muts
        LEFT JOIN
          new_muts
        USING
          (chrom,
            pos,
            ref,
            alt)
        WHERE
          new_muts.chrom IS null
    """

    destination_table = f"{dataset_id}.only_old_old"

    job_config = bigquery.QueryJobConfig(
        destination=destination_table,
        write_disposition="WRITE_APPEND",
    )

    print(f"Inserting into only_old_old for {sample_id}")
    query_job = client.query(query, job_config=job_config)
    query_job.result()
