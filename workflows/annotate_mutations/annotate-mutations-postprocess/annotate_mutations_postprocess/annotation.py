import json

import pandas as pd


def annotate_vcf(
    df: pd.DataFrame, oncogenes: set[str], tumor_suppressor_genes: set[str]
) -> pd.DataFrame:
    df_orig = df.copy()
    df = df_orig.copy()

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
