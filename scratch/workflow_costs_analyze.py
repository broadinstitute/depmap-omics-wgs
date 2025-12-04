import pandas as pd

### new
wgs = pd.read_parquet("./data/costs/depmap-omics-wgs.parquet")

wgs_cat_map = {
    "align_wgs_sample": "alignment",
    "annotate_mutations": "mutations",
    "annotate_mutations_merge": "mutations",
    "annotate_structural_variants": "SVs",
    # "bam_to_cram",
    "call_cnvs": "CNVs",
    "call_mutations": "mutations",
    "call_structural_variants": "SVs",
    "infer_msi_status": "msi",
    # "make_guide_mutation_beds": "guides",
    "prep_annotations": "mutations",
    "select_somatic_variants": "mutations",
    "select_structural_variants": "SVs",
}

wgs = wgs.merge(
    pd.DataFrame.from_dict(wgs_cat_map, orient="index", columns=["cat"])
    .reset_index()
    .rename(columns={"index": "workflow_name"}),
    how="inner",
    on="workflow_name",
)

wgs_costs = wgs.groupby(["terra_entity_name", "cat"])["cost"].agg("sum").reset_index()

wgs_costs_agg = wgs_costs.groupby("cat")["cost"].agg("sum").reset_index()
wgs_costs_n_samples = wgs_costs.value_counts("cat").reset_index()
wgs_costs_final = wgs_costs_agg.merge(wgs_costs_n_samples, how="inner", on="cat")
wgs_costs_final["cost_per_sample"] = wgs_costs_final["cost"] / wgs_costs_final["count"]

wgs_costs_final.loc[
    wgs_costs_final["cat"].isin(["CNVs", "SVs", "msi", "mutations"]), "cat"
] = "mutations+CNVs+SVs"

wgs_costs_final = (
    wgs_costs_final.groupby("cat")["cost_per_sample"].agg("sum").reset_index()
)

### old
wgs_cn = pd.read_parquet("./data/costs/DepMap_WGS_CN.parquet")

wgs_cn_cat_map = {
    "Manta_SomaticSV_v1_0": "mutations+CNVs+SVs",
    "PureCN_dev_xybAdJHsaIw": "mutations+CNVs+SVs",
    "WGS_pipeline_Z6Xikld34o8": "mutations+CNVs+SVs",
    "WGS_preprocessing": "alignment",
    "WGS_preprocessing_DRAGEN_24Q2_swqznPvpjbc": "alignment",
    "WGS_preprocessing_DRAGEN_24Q4": "alignment",
    "WGS_preprocessing_DRAGEN_24Q4_-BdjeKwv5p4": "alignment",
    "WGS_preprocessing_DRAGEN_24Q4_OJeaZWrevKw": "alignment",
    "WGS_preprocessing_WsC2mCkRJJU": "alignment",
    "WGS_preprocessing_YTIEnWr8Mqo": "alignment",
    "vcf_to_depmap_MvG9Znahp4c": "mutations+CNVs+SVs",
    "vep_sv": "mutations+CNVs+SVs",
}

wgs_cn = wgs_cn.merge(
    pd.DataFrame.from_dict(wgs_cn_cat_map, orient="index", columns=["cat"])
    .reset_index()
    .rename(columns={"index": "method_configuration_name"}),
    how="inner",
    on="method_configuration_name",
)

wgs_cn_costs = (
    wgs_cn.groupby(["terra_entity_name", "cat"])["cost"].agg("sum").reset_index()
)

wgs_cn_costs_agg = wgs_cn_costs.groupby("cat")["cost"].agg("sum").reset_index()
wgs_cn_costs_n_samples = wgs_cn_costs.value_counts("cat").reset_index()
wgs_cn_costs_final = wgs_cn_costs_agg.merge(
    wgs_cn_costs_n_samples, how="inner", on="cat"
)
wgs_cn_costs_final["cost_per_sample"] = (
    wgs_cn_costs_final["cost"] / wgs_cn_costs_final["count"]
)
