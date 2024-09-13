import json
import logging

import pandas as pd
import pysam

from annotate_mutations_postprocess.gc import calc_gc_percentage


def annotate_vcf(
    df: pd.DataFrame,
    oncogenes: set[str],
    tumor_suppressor_genes: set[str],
    fasta_path: str,
) -> pd.DataFrame:
    logging.info("Annotating VCF")

    # todo: brca1
    # transcript_likely_lof
    has_revel = df["info__oc_revel_all"].notna()

    # extract transcript IDs (field 1) given score cutoff (field 2) from values like
    # `[["ENST00000379410",0.042,0.11227],["ENST00000...`
    df.loc[has_revel, "post__transcript_likely_lof"] = (
        df.loc[has_revel, "info__oc_revel_all"]
        .apply(json.loads)
        .apply(lambda x: ";".join([y[0] for y in x if y[1] >= 0.7]))
    ).replace({"": pd.NA})

    df["post__transcript_likely_lof"] = df["post__transcript_likely_lof"].astype(
        "string"
    )

    # oncogenes and tumor suppressors
    df["post__oncogene_high_impact"] = df["info__csq__impact"].eq("HIGH") & df[
        "info__csq__symbol"
    ].isin(oncogenes)

    df["post__tumor_suppressor_high_impact"] = df["info__csq__impact"].eq("HIGH") & df[
        "info__csq__symbol"
    ].isin(tumor_suppressor_genes)

    logging.info("Calculating GC content")

    with pysam.FastaFile(fasta_path) as fasta_handle:
        df["post__gc_percentage"] = df.apply(
            lambda x: calc_gc_percentage(
                chrom=x["chromosome"],
                pos=x["position"],
                ref=x["ref"],
                variant_class=x["info__csq__variant_class"],
                window_size=200,
                fasta_handle=fasta_handle,
            ),
            axis=1,
        )

    return df
