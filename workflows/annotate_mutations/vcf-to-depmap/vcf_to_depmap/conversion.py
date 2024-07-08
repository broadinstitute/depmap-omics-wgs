from pathlib import Path

from vcf_to_depmap.utils import get_vcf_info_and_format_dtypes, read_vcf


def convert(
    vcf_path: Path,
    oncogenes_list_path: Path,
    tsg_list_path: Path,
    out_path: Path,
    force_keep: set[str],
):
    with open(oncogenes_list_path, "r") as f:
        oncogenes = set(line.strip() for line in f)

    with open(tsg_list_path, "r") as f:
        tumor_suppressor_genes = set(line.strip() for line in f)

    info_and_format_dtypes = get_vcf_info_and_format_dtypes(vcf_path)
    vcf = read_vcf(vcf_path, info_and_format_dtypes)

    pass
