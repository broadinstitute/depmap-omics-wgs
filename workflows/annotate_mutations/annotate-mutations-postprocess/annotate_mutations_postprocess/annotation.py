import json
import logging
import os
from pathlib import Path

import duckdb
import pandas as pd


def annotate_vcf(
    db_path: Path,
    parquet_dir_path: Path,
    oncogenes: set[str],
    tumor_suppressor_genes: set[str],
) -> pd.DataFrame:
    logging.info("Annotating VCF")

    try:
        os.remove(db_path)
    except OSError:
        pass

    with duckdb.connect(db_path) as db:
        logging.info(f"Reading schema and Parquet files from {parquet_dir_path}")
        db.sql(f"IMPORT DATABASE '{parquet_dir_path}'")

        # todo: brca1
        brca1 = db.sql("""
            WITH brca1_long AS (
                SELECT
                    vid,
                    k,
                    v_varchar,
                    v_float
                FROM
                    info
                WHERE
                    k IN ('oc_brca1_func_assay_class', 'oc_brca1_func_assay_score')
            ),
            assay_class AS (
                PIVOT
                    brca1_long
                ON
                    k IN ('oc_brca1_func_assay_class')
                USING
                    first(v_varchar)
                GROUP BY
                    vid
            ),
            assay_score AS (
                PIVOT
                    brca1_long
                ON
                    k IN ('oc_brca1_func_assay_score')
                USING
                    first(v_float)
                GROUP BY
                    vid
            )
            SELECT
                assay_class.vid,
                assay_class.oc_brca1_func_assay_class,
                assay_score.oc_brca1_func_assay_score
            FROM
                assay_class
            FULL OUTER JOIN
                assay_score
            ON
                assay_class.vid = assay_score.vid
        """)

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

    return df
