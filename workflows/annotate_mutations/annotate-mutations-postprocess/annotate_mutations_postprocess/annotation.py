import json
import logging
import os
from pathlib import Path

import duckdb
import pandas as pd
import pysam

from annotate_mutations_postprocess.gc import calc_gc_percentage


def annotate_vcf(
    db_path: Path,
    parquet_dir_path: Path,
    oncogenes: set[str],
    tumor_suppressor_genes: set[str],
    fasta_path: str,
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

        logging.info("Calculating GC content")
        df_gc = db.sql("""
            SELECT
                vid,
                chrom,
                pos,
                ref,
                CASE
                    WHEN length(ref) > length(alt) THEN 'deletion'
                    WHEN length(ref) > 1 AND length(alt) > 1 THEN 'substitution'
                    ELSE NULL
                END AS variant_class
            FROM
                variants
        """).df()

        with pysam.FastaFile(fasta_path) as fasta_handle:
            df_gc["gc_percentage"] = df_gc.apply(
                lambda x: calc_gc_percentage(
                    chrom=x["chrom"],
                    pos=x["pos"],
                    ref=x["ref"],
                    variant_class=x["variant_class"],
                    window_size=200,
                    fasta_handle=fasta_handle,
                ),
                axis=1,
            )

        db.sql("""
            select 
                count(distinct block_id) as num_blocks,
                count(distinct block_id) * (select block_size from pragma_database_size()) as num_bytes
                from pragma_storage_info(variants) group by all
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
