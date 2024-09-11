from pprint import pprint as pp

import pandas as pd

from annotate_mutations_postprocess.utils import cs, echo, normalize_text


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
    df["structural_relation"] = df["structural_relation"].astype("string")

    # dna repair
    dna_repair_genes = prep_dna_repair_df(dna_repair_genes)

    df = df.merge(
        dna_repair_genes, how="left", on="info__funcotation__gencode_43_hugo_symbol"
    )

    # TODO: oncokb
    df["associated_with"] = pd.NA
    df["info__oncokb_muteff"] = df["info__oncokb_muteff"].apply(normalize_text)

    df["lof"] = df["info__oncokb_muteff"].eq("loss-of-function")
    df["likely_lof"] = df["info__oncokb_muteff"].eq("likely loss-of-function")
    df["gof"] = df["info__oncokb_muteff"].eq("gain-of-function")
    df["likely_gof"] = df["info__oncokb_muteff"].eq("likely gain-of-function")

    df["driver"] = ~df["filter__multiallelic"] & df["info__civic_score"].ge(8)
    # todo: brca1

    # transcript_likely_lof
    has_revel = df["info__oc_revel_all"].notna()

    # extract transcript IDs (field 1) given score cutoff (field 2) from values like
    # `[["ENST00000379410",0.042,0.11227],["ENST00000...`
    df.loc[has_revel, "transcript_likely_lof"] = (
        df.loc[has_revel, "info__oc_revel_all"]
        .apply(json.loads)
        .apply(lambda x: ";".join([y[0] for y in x if y[1] >= 0.7]))
    ).replace({"": pd.NA})

    df["transcript_likely_lof"] = df["transcript_likely_lof"].astype("string")

    # oncogenes and tumor suppressors
    df["oncogene_high_impact"] = df["info__csq__impact"].eq("HIGH") & df[
        "info__csq__symbol"
    ].isin(oncogenes)
    df["tumor_suppressor_high_impact"] = df["info__csq__impact"].eq("HIGH") & df[
        "info__csq__symbol"
    ].isin(tumor_suppressor_genes)

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
