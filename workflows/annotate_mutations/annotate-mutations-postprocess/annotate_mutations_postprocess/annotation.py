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

        revel = db.sql("""
            WITH exploded AS (
                SELECT
                    vid,
                    unnest(json_transform_strict(v_varchar, '["json"]')) AS revel_cols
                FROM
                    info
                WHERE
                    k = 'oc_revel_all'
            ),
            split_to_cols AS (
                SELECT
                    vid,
                    json_extract_string(revel_cols, '$[0]') AS transcript_id,
                    json_extract(revel_cols, '$[1]')::DOUBLE AS score
                FROM exploded
            )
            SELECT
                vid,
                string_agg(transcript_id, ';') AS revel_transcript_ids
            FROM
                split_to_cols
            WHERE
                score >= 0.7
            GROUP BY
                vid
        """)

    # oncogenes and tumor suppressors
    df["post__oncogene_high_impact"] = df["info__csq__impact"].eq("HIGH") & df[
        "info__csq__symbol"
    ].isin(oncogenes)

    df["post__tumor_suppressor_high_impact"] = df["info__csq__impact"].eq("HIGH") & df[
        "info__csq__symbol"
    ].isin(tumor_suppressor_genes)

    return df
