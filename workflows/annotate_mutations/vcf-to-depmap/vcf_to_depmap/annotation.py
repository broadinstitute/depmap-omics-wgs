from pprint import pprint as pp

import pandas as pd

from vcf_to_depmap.utils import cs


def annotate_vcf(
    df: pd.DataFrame,
    dna_repair_genes: pd.DataFrame,
    oncogenes: set[str],
    tumor_suppressor_genes: set[str],
) -> pd.DataFrame:
    # loss of function
    df["likely_lof"] = df[cs(df, "info__csq__impact")].eq("HIGH").any(axis=1)

    # TODO: structural_relation

    dna_repair_genes = prep_dna_repair_df(dna_repair_genes)

    df = df.merge(
        dna_repair_genes, how="left", on="info__funcotation__gencode_43_hugo_symbol"
    )

    return df


def prep_dna_repair_df(dna_repair_genes: pd.DataFrame) -> pd.DataFrame:
    dna_repair_genes["gene_name"] = dna_repair_genes["gene_name"].str.split("(")
    dna_repair_genes = dna_repair_genes.explode("gene_name")
    dna_repair_genes["gene_name"] = dna_repair_genes["gene_name"].str.strip("( )")
    dna_repair_genes = dna_repair_genes.rename(
        columns={"gene_name": "info__funcotation__gencode_43_hugo_symbol"}
    )

    dna_repair_genes["dna_repair"] = (
        dna_repair_genes["accession_number"] + ": " + dna_repair_genes["activity"]
    )

    return dna_repair_genes[["info__funcotation__gencode_43_hugo_symbol", "dna_repair"]]
