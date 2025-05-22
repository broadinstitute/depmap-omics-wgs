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

query = """
    WITH
      only_old_old AS (
      SELECT
        sample_id,
        chrom,
        pos,
        ref,
        alt
      FROM
        `depmap-omics.mut_eval.only_old_old`),
      all_new AS (
      SELECT
        *
      FROM
        `depmap-omics.mut_eval.new`)
    SELECT
      only_old_old.*,
      all_new.vid,
      all_new.qual,
      all_new.ref_count,
      all_new.alt_count,
      all_new.af,
      all_new.dp,
      all_new.gt,
      all_new.ps,
      all_new.somatic,
      all_new.impactful_splice_event,
      all_new.in_clustered_event,
      all_new.segdup,
      all_new.repeat_masker,
      all_new.low_pop_prevalence,
      all_new.pon,
      all_new.rescued,
      all_new.rescued_cmc_tier,
      all_new.rescued_hess,
      all_new.rescued_met,
      all_new.rescued_oc_brca1_func_assay_score,
      all_new.rescued_oncogene_high_impact,
      all_new.rescued_oncokb_hotspot,
      all_new.rescued_oncokb_muteff,
      all_new.rescued_oncokb_oncogenic,
      all_new.rescued_tert,
      all_new.rescued_tumor_suppressor_high_impact,
      all_new.civic_desc,
      all_new.civic_id,
      all_new.civic_score,
      all_new.cosmic_tier,
      all_new.gc_prop,
      all_new.hess_driver,
      all_new.hess_signature,
      all_new.lof,
      all_new.mc,
      all_new.nmd,
      all_new.oc_brca1_func_assay_score,
      all_new.oc_gtex_gtex_gene,
      all_new.oc_gwas_catalog_disease,
      all_new.oc_gwas_catalog_pmid,
      all_new.oc_pharmgkb_id,
      all_new.oc_provean_prediction,
      all_new.oc_revel_score,
      all_new.oncogene_high_impact,
      all_new.oncokb_hotspot,
      all_new.oncokb_muteff,
      all_new.oncokb_oncogenic,
      all_new.protein_changed,
      all_new.rs,
      all_new.tumor_suppressor_high_impact,
      all_new.vep_am_class,
      all_new.vep_am_pathogenicity,
      all_new.vep_biotype,
      all_new.vep_clin_sig,
      all_new.vep_consequence,
      all_new.vep_ensp,
      all_new.vep_existing_variation,
      all_new.vep_exon,
      all_new.vep_feature,
      all_new.vep_gene,
      all_new.vep_gnom_ade_af,
      all_new.vep_gnom_adg_af,
      all_new.vep_hgnc_id,
      all_new.vep_hgvsc,
      all_new.vep_hgvsp,
      all_new.vep_impact,
      all_new.vep_intron,
      all_new.vep_loftool,
      all_new.vep_mane_select,
      all_new.vep_pli_gene_value,
      all_new.vep_poly_phen,
      all_new.vep_sift,
      all_new.vep_somatic,
      all_new.vep_swissprot,
      all_new.vep_symbol,
      all_new.vep_uniprot_isoform,
      all_new.vep_variant_class,
      all_new.vep_vid
    FROM
      only_old_old
    LEFT JOIN
      all_new
    USING
      (sample_id,
        chrom,
        pos,
        ref,
        alt)
"""

destination_table = f"{dataset_id}.only_old_new"

job_config = bigquery.QueryJobConfig(
    destination=destination_table,
    write_disposition="WRITE_APPEND",
)

print(f"Inserting into only_old_new")
query_job = client.query(query, job_config=job_config)
query_job.result()

query = """
    SELECT
      'gained' AS which,
      vid,
      sample_id,
      chrom,
      pos,
      ref,
      alt,
      qual,
      ref_count,
      alt_count,
      af,
      dp,
      gt,
      ps,
      somatic,
      impactful_splice_event,
      in_clustered_event,
      segdup,
      repeat_masker,
      low_pop_prevalence,
      pon,
      rescued,
      rescued_cmc_tier,
      rescued_hess,
      rescued_met,
      rescued_oc_brca1_func_assay_score,
      rescued_oncogene_high_impact,
      rescued_oncokb_hotspot,
      rescued_oncokb_muteff,
      rescued_oncokb_oncogenic,
      rescued_tert,
      rescued_tumor_suppressor_high_impact,
      civic_desc,
      civic_id,
      civic_score,
      cosmic_tier,
      gc_prop,
      hess_driver,
      hess_signature,
      lof,
      mc,
      nmd,
      oc_brca1_func_assay_score,
      oc_gtex_gtex_gene,
      oc_gwas_catalog_disease,
      oc_gwas_catalog_pmid,
      oc_pharmgkb_id,
      oc_provean_prediction,
      oc_revel_score,
      oncogene_high_impact,
      oncokb_hotspot,
      oncokb_muteff,
      oncokb_oncogenic,
      protein_changed,
      rs,
      tumor_suppressor_high_impact,
      vep_am_class,
      vep_am_pathogenicity,
      vep_biotype,
      vep_clin_sig,
      vep_consequence,
      vep_ensp,
      vep_existing_variation,
      vep_exon,
      vep_feature,
      vep_gene,
      vep_gnom_ade_af,
      vep_gnom_adg_af,
      vep_hgnc_id,
      vep_hgvsc,
      vep_hgvsp,
      vep_impact,
      vep_intron,
      vep_loftool,
      vep_mane_select,
      vep_pli_gene_value,
      vep_poly_phen,
      vep_sift,
      vep_somatic,
      vep_swissprot,
      vep_symbol,
      vep_uniprot_isoform,
      vep_variant_class,
      vep_vid
    FROM
      `depmap-omics.mut_eval.only_new_new`
    UNION ALL
    SELECT
      'lost' AS which,
      vid,
      sample_id,
      chrom,
      pos,
      ref,
      alt,
      qual,
      ref_count,
      alt_count,
      af,
      dp,
      gt,
      ps,
      somatic,
      impactful_splice_event,
      in_clustered_event,
      segdup,
      repeat_masker,
      low_pop_prevalence,
      pon,
      rescued,
      rescued_cmc_tier,
      rescued_hess,
      rescued_met,
      rescued_oc_brca1_func_assay_score,
      rescued_oncogene_high_impact,
      rescued_oncokb_hotspot,
      rescued_oncokb_muteff,
      rescued_oncokb_oncogenic,
      rescued_tert,
      rescued_tumor_suppressor_high_impact,
      civic_desc,
      civic_id,
      civic_score,
      cosmic_tier,
      gc_prop,
      hess_driver,
      hess_signature,
      lof,
      mc,
      nmd,
      oc_brca1_func_assay_score,
      oc_gtex_gtex_gene,
      oc_gwas_catalog_disease,
      oc_gwas_catalog_pmid,
      oc_pharmgkb_id,
      oc_provean_prediction,
      oc_revel_score,
      oncogene_high_impact,
      oncokb_hotspot,
      oncokb_muteff,
      oncokb_oncogenic,
      protein_changed,
      rs,
      tumor_suppressor_high_impact,
      vep_am_class,
      vep_am_pathogenicity,
      vep_biotype,
      vep_clin_sig,
      vep_consequence,
      vep_ensp,
      vep_existing_variation,
      vep_exon,
      vep_feature,
      vep_gene,
      vep_gnom_ade_af,
      vep_gnom_adg_af,
      vep_hgnc_id,
      vep_hgvsc,
      vep_hgvsp,
      vep_impact,
      vep_intron,
      vep_loftool,
      vep_mane_select,
      vep_pli_gene_value,
      vep_poly_phen,
      vep_sift,
      vep_somatic,
      vep_swissprot,
      vep_symbol,
      vep_uniprot_isoform,
      vep_variant_class,
      vep_vid
    FROM
      `depmap-omics.mut_eval.only_old_new`
"""

destination_table = f"{dataset_id}.diff"

job_config = bigquery.QueryJobConfig(
    destination=destination_table,
    write_disposition="WRITE_APPEND",
)

print(f"Inserting into diff")
query_job = client.query(query, job_config=job_config)
query_job.result()
