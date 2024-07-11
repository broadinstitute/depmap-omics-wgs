from pathlib import Path

from vcf_to_depmap.annotation import annotate_vcf
from vcf_to_depmap.vcf_utils import clean_vcf, get_vcf_info_and_format_dtypes, read_vcf


def convert(
    vcf_path: Path,
    oncogenes_list_path: Path,
    tsg_list_path: Path,
    out_path: Path,
    drop_cols: set[str],
    force_keep: set[str],
    compound_info_fields: set[str],
    url_encoded_col_name_regex: str | None,
    funco_sanitized_col_name_regex: str | None,
):
    with open(oncogenes_list_path, "r") as f:
        oncogenes = set(line.strip() for line in f)

    with open(tsg_list_path, "r") as f:
        tumor_suppressor_genes = set(line.strip() for line in f)

    info_and_format_dtypes = get_vcf_info_and_format_dtypes(
        vcf_path, compound_info_fields
    )

    df = read_vcf(vcf_path)

    df = clean_vcf(
        df,
        info_and_format_dtypes,
        drop_cols,
        url_encoded_col_name_regex,
        funco_sanitized_col_name_regex,
    )

    annotate_vcf(df, oncogenes, tumor_suppressor_genes)

    pass
