from pprint import pprint as pp

import pandas as pd

from annotate_mutations_postprocess.utils import cs


def annotate_vcf(
    df: pd.DataFrame,
    dna_repair_genes: pd.DataFrame,
    oncogenes: set[str],
    tumor_suppressor_genes: set[str],
    cosmic_fusions: pd.DataFrame,
    cosmic_translocation_partners: pd.DataFrame,
) -> pd.DataFrame:
    df_orig = df.copy()
    df = df_orig.copy()

    # loss of function
    df["likely_lof"] = df[cs(df, "info__csq__impact")].eq("HIGH").any(axis=1)

    # structural relation
    df = (
        df.merge(
            cosmic_fusions,
            how="left",
            left_on="info__funcotation__gencode_43_hugo_symbol",
            right_on="gene",
        )
        .merge(cosmic_translocation_partners, how="left", on="gene")
        .drop(columns="gene")
    )

    has_structural_relation = (
        df["translocation_partners"].isna() & df["fusion_genes"].notna()
    )

    df.loc[has_structural_relation, "structural_relation"] = df.loc[
        has_structural_relation, "fusion_genes"
    ].str.extract(r"^.+::.+\(([A-Z 0-9]+)\):", expand=False)

    dna_repair_genes = prep_dna_repair_df(dna_repair_genes)

    df = df.merge(
        dna_repair_genes, how="left", on="info__funcotation__gencode_43_hugo_symbol"
    )

    # TODO: oncokb

    df["driver"] = ~df["filter__multiallelic"] & df["info__civic_score"].ge(8)

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
