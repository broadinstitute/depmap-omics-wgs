import json
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
from google.cloud import storage
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

sample_sets = ws.get_entities("sample_set")

sample_ids = [
    x["entityName"]
    for x in sample_sets.loc[
        sample_sets["sample_set_id"].eq("warp_eval"), "samples"
    ].values[0]["items"]
]

samples = ws.get_entities("sample")

eval_samples = samples.loc[
    samples["sample_id"].isin(sample_ids),
    [
        "sample_id",
        "mut_annot_vcf",
        "mut_annot_vcf_index",
        "mut_duckdb",
        "mut_somatic_variants",
        "wgs_cn_mutect2_filtered_vcf",
        "wgs_cn_mutect2_filtered_vcf_idx",
        "wgs_cn_full_annotated_parquet",
        "wgs_cn_depmap_maf_25q2",
    ],
]

eval_samples.set_index("sample_id", inplace=True)

eval_samples = eval_samples.replace({"": pd.NA})
eval_samples.loc["CDS-00Nrci", "mut_somatic_variants"] = (
    "gs://fc-secure-ee51f36b-a39f-4191-93a2-940ae7960793/submissions/final-outputs/6e31bc09-1af7-4a0b-b755-28aec146c0fd/postprocess_mutations/8fe872bc-0b37-45e3-a346-1e13876acf11/call-postprocess/attempt-2/somatic_variants.parquet"
)
eval_samples.loc["CDS-00Nrci", "mut_duckdb"] = """{
    "itemsType": "AttributeValue",
    "items": [
        "gs://fc-secure-ee51f36b-a39f-4191-93a2-940ae7960793/submissions/final-outputs/cf6d7a90-e64e-45a1-a105-154cc1a9b31a/annotate_mutations_merge/8618abee-0705-4c9c-9a14-7716e2aa2f4e/call-vcf_to_duckdb/glob-f173131bbcb257c82200c4ab7733d446/load.sql",
        "gs://fc-secure-ee51f36b-a39f-4191-93a2-940ae7960793/submissions/final-outputs/cf6d7a90-e64e-45a1-a105-154cc1a9b31a/annotate_mutations_merge/8618abee-0705-4c9c-9a14-7716e2aa2f4e/call-vcf_to_duckdb/glob-f173131bbcb257c82200c4ab7733d446/schema.sql",
        "gs://fc-secure-ee51f36b-a39f-4191-93a2-940ae7960793/submissions/final-outputs/cf6d7a90-e64e-45a1-a105-154cc1a9b31a/annotate_mutations_merge/8618abee-0705-4c9c-9a14-7716e2aa2f4e/call-vcf_to_duckdb/glob-f173131bbcb257c82200c4ab7733d446/val_info_types.parquet",
        "gs://fc-secure-ee51f36b-a39f-4191-93a2-940ae7960793/submissions/final-outputs/cf6d7a90-e64e-45a1-a105-154cc1a9b31a/annotate_mutations_merge/8618abee-0705-4c9c-9a14-7716e2aa2f4e/call-vcf_to_duckdb/glob-f173131bbcb257c82200c4ab7733d446/vals_info.parquet",
        "gs://fc-secure-ee51f36b-a39f-4191-93a2-940ae7960793/submissions/final-outputs/cf6d7a90-e64e-45a1-a105-154cc1a9b31a/annotate_mutations_merge/8618abee-0705-4c9c-9a14-7716e2aa2f4e/call-vcf_to_duckdb/glob-f173131bbcb257c82200c4ab7733d446/variants.parquet"
    ]
}"""
eval_samples["mut_duckdb"] = eval_samples["mut_duckdb"].apply(
    lambda x: json.loads(x) if x is not pd.NA else pd.NA
)


def cache_gcs_blob(gcs_url: str, local_dir: Path, client: storage.Client) -> None:
    gcs_url = gcs_url.strip()
    assert gcs_url.startswith("gs://")

    bucket_name, *blob_parts = gcs_url.replace("gs://", "").split("/")
    blob_path = "/".join(blob_parts)
    filename = Path(blob_path).name
    local_file = local_dir / filename

    if not local_file.exists():
        bucket = client.bucket(bucket_name)
        blob = bucket.blob(blob_path)
        blob.download_to_filename(str(local_file))


def cache_output_files(c: str, s: pd.Series) -> None:
    local_dir = Path(f"./data/mut_eval/{c}")
    local_dir.mkdir(parents=True, exist_ok=True)

    client = storage.Client(project="depmap-omics")

    with ThreadPoolExecutor(max_workers=10) as executor:
        for sample_id, x in s.items():
            sample_dir = local_dir / str(sample_id)
            sample_dir.mkdir(parents=True, exist_ok=True)

            if isinstance(x, dict):
                for y in x["items"]:
                    executor.submit(
                        cache_gcs_blob, gcs_url=y, local_dir=sample_dir, client=client
                    )
            else:
                executor.submit(
                    cache_gcs_blob, gcs_url=x, local_dir=sample_dir, client=client
                )

        executor.shutdown(wait=True)


for col in [
    # "mut_annot_vcf",
    # "mut_annot_vcf_index",
    # "mut_duckdb",
    "mut_somatic_variants",
    # "wgs_cn_mutect2_filtered_vcf",
    # "wgs_cn_mutect2_filtered_vcf_idx",
    # "wgs_cn_full_annotated_parquet",
    "wgs_cn_depmap_maf_25q2",
]:
    print(col)
    cache_output_files(c=col, s=eval_samples[col])


def load_somatic_variants_old() -> pd.DataFrame:
    cache_path = Path("./data/mut_eval/somatic_variants_old.parquet")

    if cache_path.exists():
        return pd.read_parquet(cache_path)

    local_dir = Path("./data/mut_eval/wgs_cn_depmap_maf_25q2")

    def download_and_read(f: Path, sample_id: str) -> pd.DataFrame:
        return pd.read_csv(f, dtype="string").assign(sample_id=sample_id)

    somatic_variants_old_arr = []

    with ThreadPoolExecutor(max_workers=10) as executor:
        futures = []

        for d in local_dir.iterdir():
            if d.is_dir():
                sample_id = d.name

                for csvgz in d.glob("*.csv.gz"):
                    futures.append(
                        executor.submit(download_and_read, f=csvgz, sample_id=sample_id)
                    )

        for future in as_completed(futures):
            df = future.result()
            somatic_variants_old_arr.append(df)

    df = pd.concat(somatic_variants_old_arr, ignore_index=True)

    df = df.drop(columns=["031225"])

    bool_cols = [
        "oncogene_high_impact",
        "tumor_suppressor_high_impact",
        "lof",
        "driver",
        "likely_driver",
        "likely_gof",
        "likely_lof",
        "hess_driver",
        "segdup",
        "rescue",
    ]

    df[bool_cols] = df[bool_cols].astype("string")
    df[bool_cols] = df[bool_cols].isin(["True", "Y"]).astype("boolean")

    df = df.astype(
        {
            "chrom": "string",
            "pos": "int64",
            "ref": "string",
            "alt": "string",
            "af": "float64",
            "dp": "int64",
            "ref_count": "int64",
            "alt_count": "int64",
            "gt": "string",
            "ps": "Int64",
            "variant_type": "string",
            "variant_info": "string",
            "dna_change": "string",
            "protein_change": "string",
            "hugo_symbol": "string",
            "exon": "string",
            "intron": "string",
            "ensembl_gene_id": "string",
            "ensembl_feature_id": "string",
            "hgnc_name": "string",
            "hgnc_family": "string",
            "uniprot_id": "string",
            "dbsnp_rs_id": "string",
            "gc_content": "float64",
            "lof_gene_name": "string",
            "lof_gene_id": "string",
            "lof_number_of_transcripts_in_gene": "Int64",
            "lof_percent_of_transcripts_affected": "Float64",
            "nmd": "string",
            "clnsig": "string",
            "molecular_consequence": "string",
            "af_exac": "Float64",
            "af_esp": "Float64",
            "af_tgp": "Float64",
            "vep_impact": "string",
            "vep_biotype": "string",
            "vep_hgnc_id": "string",
            "vep_existing_variation": "string",
            "vep_mane_select": "string",
            "vep_ensp": "string",
            "vep_swissprot": "string",
            "sift": "string",
            "polyphen": "string",
            "gnomade_af": "Float64",
            "gnomadg_af": "Float64",
            "am_class": "string",
            "am_pathogenicity": "Float64",
            "vep_clin_sig": "string",
            "vep_somatic": "string",
            "vep_pli_gene_value": "Float64",
            "vep_loftool": "Float64",
            "achilles_top_genes": "string",
            "structural_relation": "string",
            "associated_with": "string",
            "transcript_likely_lof": "string",
            "brca1_func_score": "Float64",
            "civic_id": "Float64",
            "civic_description": "string",
            "civic_score": "Float64",
            "popaf": "Float64",
            "hess_signture": "string",
            "revel_score": "Float64",
            "pharmgkb_id": "string",
            "dida_id": "string",
            "dida_name": "string",
            "gwas_disease": "string",
            "gwas_pmid": "Int64",
            "gtex_gene": "string",
            "cosmic_tier": "Int64",
            "oncokb_effect": "string",
            "oncokb_hotspot": "string",
            "oncokb_oncogenic": "string",
            "provean_prediction": "string",
            "spliceai_ds_ag": "Float64",
            "spliceai_ds_al": "Float64",
            "spliceai_ds_dg": "Float64",
            "spliceai_ds_dl": "Float64",
            "spliceai_dp_ag": "Float64",
            "spliceai_dp_al": "Float64",
            "spliceai_dp_dg": "Float64",
            "spliceai_dp_dl": "Float64",
            "sample_id": "string",
            "oncogene_high_impact": "boolean",
            "tumor_suppressor_high_impact": "boolean",
            "lof": "boolean",
            "driver": "boolean",
            "likely_driver": "boolean",
            "likely_gof": "boolean",
            "likely_lof": "boolean",
            "hess_driver": "boolean",
            "segdup": "boolean",
            "rm": "boolean",
            "rescue": "boolean",
        }
    )

    df["civic_id"] = df["civic_id"].astype("Int64")

    df = df.loc[
        :,
        [
            "sample_id",
            "chrom",
            "pos",
            "ref",
            "alt",
            "af",
            "dp",
            "ref_count",
            "alt_count",
            "gt",
            "ps",
            "variant_type",
            "variant_info",
            "dna_change",
            "protein_change",
            "hugo_symbol",
            "exon",
            "intron",
            "ensembl_gene_id",
            "ensembl_feature_id",
            "hgnc_name",
            "hgnc_family",
            "uniprot_id",
            "dbsnp_rs_id",
            "achilles_top_genes",
            "af_esp",
            "af_exac",
            "af_tgp",
            "am_class",
            "am_pathogenicity",
            "associated_with",
            "brca1_func_score",
            "civic_description",
            "civic_id",
            "civic_score",
            "clnsig",
            "cosmic_tier",
            "dida_id",
            "dida_name",
            "driver",
            "gc_content",
            "gnomade_af",
            "gnomadg_af",
            "gtex_gene",
            "gwas_disease",
            "gwas_pmid",
            "hess_driver",
            "hess_signture",
            "likely_driver",
            "likely_gof",
            "likely_lof",
            "lof",
            "lof_gene_id",
            "lof_gene_name",
            "lof_number_of_transcripts_in_gene",
            "lof_percent_of_transcripts_affected",
            "molecular_consequence",
            "nmd",
            "oncogene_high_impact",
            "oncokb_effect",
            "oncokb_hotspot",
            "oncokb_oncogenic",
            "pharmgkb_id",
            "polyphen",
            "popaf",
            "provean_prediction",
            "rescue",
            "revel_score",
            "rm",
            "segdup",
            "sift",
            "spliceai_dp_ag",
            "spliceai_dp_al",
            "spliceai_dp_dg",
            "spliceai_dp_dl",
            "spliceai_ds_ag",
            "spliceai_ds_al",
            "spliceai_ds_dg",
            "spliceai_ds_dl",
            "structural_relation",
            "transcript_likely_lof",
            "tumor_suppressor_high_impact",
            "vep_biotype",
            "vep_clin_sig",
            "vep_ensp",
            "vep_existing_variation",
            "vep_hgnc_id",
            "vep_impact",
            "vep_loftool",
            "vep_mane_select",
            "vep_pli_gene_value",
            "vep_somatic",
            "vep_swissprot",
        ],
    ]

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(cache_path)

    return df


def load_somatic_variants_new() -> pd.DataFrame:
    cache_path = Path("./data/mut_eval/somatic_variants_new.parquet")

    if cache_path.exists():
        return pd.read_parquet(cache_path)

    local_dir = Path("./data/mut_eval/mut_somatic_variants")

    def download_and_read(f: Path, sample_id: str) -> pd.DataFrame:
        return pd.read_parquet(f).assign(sample_id=sample_id)

    somatic_variants_new_arr = []

    with ThreadPoolExecutor(max_workers=10) as executor:
        futures = []

        for d in local_dir.iterdir():
            if d.is_dir():
                sample_id = d.name

                for parq in d.glob("*.parquet"):
                    futures.append(
                        executor.submit(download_and_read, f=parq, sample_id=sample_id)
                    )

        for future in as_completed(futures):
            somatic_variants_new_arr.append(future.result())

    df = pd.concat(somatic_variants_new_arr, ignore_index=True)

    df = df.loc[
        :,
        [
            "sample_id",
            "chrom",
            "pos",
            "ref",
            "alt",
            "af",
            "dp",
            "ref_count",
            "alt_count",
            "gt",
            "ps",
            "variant_type",
            "variant_info",
            "dna_change",
            "protein_change",
            "hugo_symbol",
            "exon",
            "intron",
            "ensembl_gene_id",
            "ensembl_feature_id",
            "hgnc_name",
            "hgnc_family",
            "uniprot_id",
            "dbsnp_rs_id",
            "am_class",
            "am_pathogenicity",
            "brca1_func_score",
            "civic_description",
            "civic_id",
            "civic_score",
            "cosmic_tier",
            "gc_content",
            "gnomade_af",
            "gnomadg_af",
            "gtex_gene",
            "gwas_disease",
            "gwas_pmid",
            "hess_driver",
            "hess_signature",
            "lof_gene_id",
            "lof_gene_name",
            "lof_number_of_transcripts_in_gene",
            "lof_prop_of_transcripts_affected",
            "molecular_consequence",
            "nmd",
            "oncogene_high_impact",
            "oncokb_effect",
            "oncokb_hotspot",
            "oncokb_oncogenic",
            "pharmgkb_id",
            "polyphen",
            "provean_prediction",
            "rescue",
            "revel_score",
            "sift",
            "transcript_likely_lof",
            "tumor_suppressor_high_impact",
            "vep_biotype",
            "vep_clin_sig",
            "vep_ensp",
            "vep_existing_variation",
            "vep_hgnc_id",
            "vep_impact",
            "vep_loftool",
            "vep_mane_select",
            "vep_pli_gene_value",
            "vep_somatic",
            "vep_swissprot",
        ],
    ]

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(cache_path)

    return df


somatic_variants_old = load_somatic_variants_old()
somatic_variants_new = load_somatic_variants_new()

somatic_variants_old = somatic_variants_old.rename(
    columns={
        "lof_percent_of_transcripts_affected": "lof_prop_of_transcripts_affected",
        "hess_signture": "hess_signature",
    }
)

# todo: temp
somatic_variants_old = somatic_variants_old.loc[
    somatic_variants_old["sample_id"].isin(somatic_variants_new["sample_id"])
]
somatic_variants_new = somatic_variants_new.loc[
    somatic_variants_new["sample_id"].isin(somatic_variants_old["sample_id"])
]

variant_cols = ["sample_id", "chrom", "pos", "ref", "alt"]

variant_comp = (
    somatic_variants_old[variant_cols]
    .merge(
        somatic_variants_new[variant_cols],
        how="outer",
        indicator=True,
    )
    .set_index(variant_cols)
)["_merge"]

only_old = somatic_variants_old.set_index(variant_cols).loc[
    variant_comp.loc[variant_comp.eq("left_only")].index
]
only_old["which"] = "old"

only_new = somatic_variants_new.set_index(variant_cols).loc[
    variant_comp.loc[variant_comp.eq("right_only")].index
]
only_new["which"] = "new"

diff_df = pd.concat([only_old.reset_index(), only_new.reset_index()], ignore_index=True)
