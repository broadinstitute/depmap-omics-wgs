import logging
from pathlib import Path

import pandas as pd

from compute_molecular_signatures.types import Maf, MutSigVariants, TypedDataFrame
from compute_molecular_signatures.utils import type_data_frame


def make_maf(muts_path: Path) -> TypedDataFrame[Maf]:
    """
    Read the mutations Parquet file, apply filters, and convert to a MAF-like data
    frame.

    :param muts_path: path to the mutations Parquet file
    :return: a MAF-like data frame
    """

    logging.info(f"Reading mutations from {muts_path}")
    df = pd.read_parquet(
        muts_path,
        columns=[
            "sample_id",
            "chrom",
            "pos",
            "ref",
            "alt",
            "alt_count",
            "af",
            "vep_gnom_ade_af",
            "vep_gnom_adg_af",
        ],
    )

    df = type_data_frame(df, MutSigVariants)

    # filter using gnomad stats (optionally imputed with 0.0 prevalence)
    logging.info("Filtering mutations")
    df = df.loc[
        df["alt_count"].ge(5)
        & df["vep_gnom_ade_af"].fillna(0.0).lt(1e-6)
        & df["vep_gnom_adg_af"].fillna(0.0).lt(1e-6)
    ].copy()

    logging.info("Making MAF data frame")
    df["Chromosome"] = df["chrom"].str.replace("^chr", "", regex=True)

    logging.info("Adjusting pos/ref/alt for insertions/deletions")
    is_ins = df["ref"].str.len().gt(df["alt"].str.len())
    is_del = df["ref"].str.len().lt(df["alt"].str.len())

    # confirm ref/alt formats are as expected
    assert bool(
        df.loc[is_ins, "ref"].str.slice(stop=1).eq(df.loc[is_ins, "alt"]).all()
    ) and bool(df.loc[is_ins, "alt"].str.len().eq(1).all()), (
        "Unexpected ref/alt formats found for insertions"
    )

    assert bool(
        df.loc[is_del, "alt"].str.slice(stop=1).eq(df.loc[is_del, "ref"]).all()
    ) and bool(df.loc[is_del, "ref"].str.len().eq(1).all()), (
        "Unexpected ref/alt formats found for deletions"
    )

    # strip shared leftmost allele
    df.loc[is_ins, "pos"] += 1
    df.loc[is_ins, "ref"] = df.loc[is_ins, "ref"].str.slice(start=1)
    df.loc[is_ins, "alt"] = "-"

    df.loc[is_del, "pos"] += 1
    df.loc[is_del, "ref"] = "-"
    df.loc[is_del, "alt"] = df.loc[is_del, "alt"].str.slice(start=1)

    # make minimal MAF-like data frame
    df = df.rename(
        columns={
            "sample_id": "Tumor_Sample_Barcode",
            "Chromosome": "Chromosome",
            "pos": "Start_Position",
            "ref": "Reference_Allele",
            "alt": "Tumor_Seq_Allele2",
        }
    )

    df["Hugo_Symbol"] = "Unknown"  # not needed for signatureanalyzer

    return type_data_frame(df, Maf, remove_unknown_cols=True)
