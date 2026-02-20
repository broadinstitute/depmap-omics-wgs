from pathlib import Path
from typing import Type

import pandas as pd

from compute_molecular_signatures.types import Maf, PanderaBaseSchema, TypedDataFrame


def make_maf(muts_path: Path) -> pd.DataFrame:
    """
    Read the mutations Parquet file, apply filters, and convert to a MAF-like data
    frame.

    :param muts_path: path to the mutations Parquet file
    :return: a MAF-like data frame
    """

    df = pd.read_parquet(muts_path)

    # filter using gnomad stats (optionally imputed with 0.0 prevalence)
    df = df.loc[
        df["alt_count"].ge(5)
        & df["vep_gnom_ade_af"].fillna(0.0).lt(1e-6)
        & df["vep_gnom_adg_af"].fillna(0.0).lt(1e-6)
    ].copy()

    df["Chromosome"] = df["chrom"].str.replace("^chr", "", regex=True)

    def normalize_indel(row: pd.Series) -> pd.Series:
        """
        Normalize (left-align) indelx and adjust positions.

        :param row: a single variant
        :return: series of (adjusted start position, ref allele, alt allele)
        """

        ref = row.ref
        alt = row.alt
        pos = row.pos

        # deletion
        if len(ref) > len(alt):
            idx = ref.find(alt)
            if idx != -1:
                return pd.Series([pos + idx + 1, ref[idx + len(alt) :], "-"])

        # insertion
        if len(ref) < len(alt):
            idx = alt.find(ref)
            if idx != -1:
                return pd.Series([pos + idx + 1, "-", alt[idx + len(ref) :]])

        # SNV or complex
        return pd.Series([pos, ref, alt])

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


def type_data_frame(
    df: pd.DataFrame,
    pandera_schema: Type[PanderaBaseSchema],
    reorder_cols: bool = True,
    remove_unknown_cols: bool = False,
) -> TypedDataFrame[PanderaBaseSchema]:
    """
    Coerce a data frame into one specified by a Pandera schema and optionally reorder
    columns and remove unknown columns.

    :param df: a data frame
    :param pandera_schema: a Pandera schema
    :param reorder_cols: reorder columns as specified in the schema
    :param remove_unknown_cols: remove columns not specified in the schema
    :return: a data frame validated with the provided Pandera schema
    """

    if len(df) == 0:
        # make an empty data frame that conforms to the Pandera schema
        s = pandera_schema.to_schema()

        # `example` doesn't know how to instantiate columns with structured data
        dict_cols = []
        list_cols = []

        for c in s.columns:
            if s.columns[c].dtype.type is dict:
                dict_cols.append(c)
                s = s.remove_columns([c])
            elif s.columns[c].dtype.type is list:
                list_cols.append(c)
                s = s.remove_columns([c])

        df = pd.DataFrame(s.example(size=1))

        if len(dict_cols) > 0:
            for c in dict_cols:
                df[c] = [{}] * len(df)

        if len(list_cols) > 0:
            for c in list_cols:
                df[c] = [[]] * len(df)

        df = df.iloc[:0]
        return TypedDataFrame[pandera_schema](df)

    if not remove_unknown_cols and not reorder_cols:
        # can type and return
        return TypedDataFrame[pandera_schema](df)

    # we need to collect the current columns and schema columns (in original orders)
    current_cols = list(df.columns)
    schema_cols = list(pandera_schema.to_schema().columns.keys())

    if remove_unknown_cols:
        # drop excess columns (if any)
        excess_cols = list(set(current_cols) - set(schema_cols))

        if len(excess_cols) > 0:
            df = df.drop(columns=excess_cols)
            current_cols = list(df.columns)

    # `df` might contain extra columns, but we can still type it now
    df = TypedDataFrame[pandera_schema](df)

    if reorder_cols:
        # put columns in schema order, with extra columns in original order at the end
        all_cols = schema_cols.copy()
        all_cols.extend(current_cols)
        all_cols = list(dict.fromkeys(all_cols))
        df = TypedDataFrame[pandera_schema](df.loc[:, all_cols])

    return df
