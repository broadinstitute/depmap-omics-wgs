import re
from pathlib import Path

import pandas as pd

from annotate_mutations_postprocess.annotation import annotate_vcf
from annotate_mutations_postprocess.vcf_utils import (
    clean_vcf,
    get_vcf_info_and_format_dtypes,
    read_vcf,
)


def convert(
    vcf_path: Path,
    dna_repair_genes_path: Path,
    oncogenes_path: Path,
    tumor_suppressor_genes_path: Path,
    cosmic_fusions_path: Path,
    cosmic_translocation_partners_path: Path,
    out_path: Path,
    drop_cols: set[str],
    na_cols: set[str],
    bool_cols: set[str],
    force_keep: set[str],
    compound_info_fields: set[str],
    url_encoded_col_name_regex: re.Pattern,
    funco_sanitized_col_name_regex: re.Pattern,
):
    dna_repair_genes = pd.read_csv(dna_repair_genes_path, dtype="string")

    with open(oncogenes_path, "r") as f:
        oncogenes = set(line.strip() for line in f)

    with open(tumor_suppressor_genes_path, "r") as f:
        tumor_suppressor_genes = set(line.strip() for line in f)

    cosmic_fusions = pd.read_csv(cosmic_fusions_path, dtype="string")
    cosmic_translocation_partners = pd.read_csv(
        cosmic_translocation_partners_path, dtype="string"
    )

    info_and_format_dtypes = get_vcf_info_and_format_dtypes(
        vcf_path, compound_info_fields
    )

    df = read_vcf(vcf_path)

    df = clean_vcf(
        df=df,
        info_and_format_dtypes=info_and_format_dtypes,
        drop_cols=drop_cols,
        na_cols=na_cols,
        bool_cols=bool_cols,
        url_encoded_col_name_regex=url_encoded_col_name_regex,
        funco_sanitized_col_name_regex=funco_sanitized_col_name_regex,
    )

    df = annotate_vcf(
        df=df,
        dna_repair_genes=dna_repair_genes,
        oncogenes=oncogenes,
        tumor_suppressor_genes=tumor_suppressor_genes,
        cosmic_fusions=cosmic_fusions,
        cosmic_translocation_partners=cosmic_translocation_partners,
    )

    pass
