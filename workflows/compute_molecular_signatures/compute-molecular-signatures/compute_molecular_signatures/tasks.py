import logging
from pathlib import Path

import pandas as pd

from compute_molecular_signatures.types import Maf, TypedDataFrame
from compute_molecular_signatures.utils import normalize_indel, type_data_frame


def make_maf(muts_path: Path) -> TypedDataFrame[Maf]:
    """
    Read the mutations Parquet file, apply filters, and convert to a MAF-like data
    frame.

    :param muts_path: path to the mutations Parquet file
    :return: a MAF-like data frame
    """

    logging.info(f"Reading mutations from {muts_path}")
    df = pd.read_parquet(muts_path)

    # filter using gnomad stats (optionally imputed with 0.0 prevalence)
    logging.info("Filtering mutations")
    df = df.loc[
        df["alt_count"].ge(5)
        & df["vep_gnom_ade_af"].fillna(0.0).lt(1e-6)
        & df["vep_gnom_adg_af"].fillna(0.0).lt(1e-6)
    ].copy()

    logging.info("Making MAF data frame")
    df["Chromosome"] = df["chrom"].str.replace("^chr", "", regex=True)

    df[["Start_Position", "Reference_Allele", "Tumor_Seq_Allele2"]] = df.apply(
        normalize_indel, axis=1
    )

    # make minimal MAF-like data frame
    df = df.rename(
        columns={
            "sample_id": "Tumor_Sample_Barcode",
            "Chromosome": "Chromosome",
            "Start_Position": "Start_Position",
            "Reference_Allele": "Reference_Allele",
            "Tumor_Seq_Allele2": "Tumor_Seq_Allele2",
        }
    )

    df["Hugo_Symbol"] = "Unknown"  # not needed for signatureanalyzer

    return type_data_frame(df, Maf, remove_unknown_cols=True)
