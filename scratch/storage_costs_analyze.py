from pprint import pp

from nebelung.terra_workspace import TerraWorkspace
from pd_flatten import pd_flatten

from depmap_omics_wgs.gcs import get_objects_metadata

### new
ws_new = TerraWorkspace("broad-firecloud-ccle", "depmap-omics-wgs")
samples_new = ws_new.get_entities("sample")

samples_new_done = (
    samples_new.loc[
        :,
        [
            "sample_id",
            # "aligned_sequencing_alignment_id",
            # "analysis_ready_bai",
            # "analysis_ready_bam",
            # "analysis_ready_crai",
            # "analysis_ready_cram",
            # "analysis_ready_cram_size",
            # "automation_status",
            # "cell_line_name",
            "cnv_bin_coverage",
            "cnv_cn_by_gene_weighted_mean",
            "cnv_input_params",
            "cnv_read_cov_bin",
            "cnv_segments",
            # "delivery_crai_bai",
            # "delivery_cram_bam",
            # "delivery_cram_bam_size",
            # "delivery_file_format",
            # "delivery_ref",
            # "delivery_ref_alt",
            # "delivery_ref_amb",
            # "delivery_ref_ann",
            # "delivery_ref_bwt",
            # "delivery_ref_dict",
            # "delivery_ref_fasta",
            # "delivery_ref_fasta_index",
            # "delivery_ref_pac",
            # "delivery_ref_sa",
            # "delivery_sequencing_alignment_id",
            "guide_bed_avana",
            "guide_bed_brunello",
            "guide_bed_humagne",
            "guide_bed_ky",
            "guide_bed_tko",
            # "model_condition_id",
            # "model_id",
            "msisensor2_output",
            "msisensor2_output_dis",
            "msisensor2_output_somatic",
            # "msisensor2_score",
            "mut_annot_bcftools_fixed_vcf",
            "mut_annot_bcftools_vcf",
            "mut_annot_open_cravat_vcf",
            "mut_annot_snpeff_snpsift_vcf",
            "mut_annot_vcf",
            "mut_annot_vcf_index",
            "mut_annot_vep_vcf",
            "mut_duckdb",
            # "mut_enriched_variants",
            "mut_filtering_stats",
            "mut_sig_variants",
            "mut_somatic_variants",
            "mut_vcf",
            "mut_vcf_idx",
            # "non_cancerous",
            # "omics_profile_id",
            # "peddep",
            "preprocess_md_metrics",
            # "ref",
            # "ref_alt",
            # "ref_amb",
            # "ref_ann",
            # "ref_bwt",
            # "ref_dict",
            # "ref_fasta",
            # "ref_fasta_index",
            # "ref_pac",
            # "ref_sa",
            # "stripped_cell_line_name",
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
        ],
    ]
    .dropna()
    .sample(n=100)
)

urls_new = (
    pd_flatten(samples_new_done)
    .drop(columns=["mut_duckdb__itemsType"])
    .melt(id_vars="sample_id")
    .drop_duplicates(subset="value")
)

urls_new = urls_new.loc[
    urls_new["value"].notna()
    & urls_new["value"].str.startswith("gs://fc-secure-ee51f36b")
]

pp(urls_new["variable"].value_counts().to_dict())

blobs_new = get_objects_metadata(urls_new["value"])
total_size_new = blobs_new["size"].sum()
size_per_sample_new = total_size_new / 100

### old
ws_old = TerraWorkspace("broad-firecloud-ccle", "DepMap_WGS_CN")
samples_old = ws_old.get_entities("sample")

samples_old_done = samples_old.sample(n=100)

urls_old = (
    pd_flatten(samples_old_done)
    .melt(id_vars="sample_id")
    .drop_duplicates(subset="value")
)

urls_old = urls_old.loc[
    urls_old["value"].notna()
    & urls_old["value"].str.startswith("gs://fc-secure-bd7b8bc9")
]

pp(urls_old["variable"].value_counts().to_dict())

blobs_old = get_objects_metadata(urls_old["value"])
total_size_old = blobs_old["size"].sum()
size_per_sample_old = total_size_old / 100
