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
    min_af: float = 0.15,
    min_depth: int = 2,
) -> None:
    logging.info("Annotating VCF")

    try:
        os.remove(db_path)
    except OSError:
        pass

    with duckdb.connect(db_path) as db:
        logging.info(f"Reading schema and Parquet files from {parquet_dir_path}")
        db.sql(f"IMPORT DATABASE '{parquet_dir_path}'")

        db.register(
            "oncogenes",
            pd.DataFrame({"gene_name": list(oncogenes)}, dtype="string"),
        )

        db.register(
            "tumor_suppressor_genes",
            pd.DataFrame({"gene_name": list(tumor_suppressor_genes)}, dtype="string"),
        )

        db.sql("""
            CREATE OR REPLACE VIEW brca1 AS (
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
            )
        """)

        db.sql("""
            CREATE OR REPLACE VIEW revel AS (
                WITH exploded AS (
                    SELECT
                        vid,
                        unnest(
                            json_transform_strict(v_varchar, '["json"]')
                        ) AS revel_cols
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
            )
        """)

        db.sql("""
            CREATE OR REPLACE VIEW vep AS (
                WITH vep_exploded AS (
                    SELECT
                        vid,
                        json_transform(
                            UNNEST(v_json_arr),
                            '{
                                "symbol": "VARCHAR",
                                "consequence": "VARCHAR",
                                "biotype": "VARCHAR",
                                "clin_sig": "VARCHAR",
                                "ensp": "VARCHAR",
                                "existing_variation": "VARCHAR",
                                "hgnc_id": "VARCHAR",
                                "impact": "VARCHAR",
                                "loftool": "DOUBLE",
                                "mane_select": "VARCHAR",
                                "pli_gene_value": "DOUBLE",
                                "somatic": "VARCHAR",
                                "swissprot": "VARCHAR"
                            }'
                        ) AS csq
                    FROM
                        info
                    WHERE
                        k = 'csq'
                )
                SELECT
                    vid,
                    csq.*
                FROM
                    vep_exploded
            )
        """)

        db.sql(f"""
            CREATE OR REPLACE VIEW vals_wide AS (
                SELECT
                    vid,
                    t_af.af,
                    t_dp.dp,
                    t_gt.gt,
                    t_ps.ps
                FROM
                    vals
                
                NATURAL JOIN (
                    SELECT
                        vid,
                        v_float AS af,
                    FROM
                        vals
                    WHERE
                        k = 'af'
                ) t_af
                
                NATURAL JOIN (
                    SELECT
                        vid,
                        v_integer AS dp
                    FROM
                        vals
                    WHERE
                        k = 'dp'
                ) t_dp
                
                NATURAL JOIN (
                    SELECT
                        vid,
                        v_varchar AS gt
                    FROM
                        vals
                    WHERE
                        k = 'gt'
                ) t_gt
                
                NATURAL JOIN (
                    SELECT
                        vid,
                        v_integer AS ps
                    FROM
                        vals
                    WHERE
                        k = 'ps'
                ) t_ps
            )
        """)

        db.sql("""
            CREATE OR REPLACE VIEW filters AS (
                WITH val_filters AS (
                    SELECT
                        vid,
                        unnest(filters) AS filter
                    FROM
                        variants
                ),
                as_filters AS (
                    SELECT
                        vid,
                        unnest(str_split_regex(v_varchar, '\\s*\\|\\s*')) AS filter
                    FROM
                        info
                    WHERE
                        k = 'as_filter_status'
                )
                SELECT
                    vid,
                    lower(filter) AS filter
                FROM
                    val_filters
                UNION 
                SELECT
                    vid,
                    lower(filter) AS filter
                FROM
                    as_filters
            )
        """)

        db.sql(f"""
            SELECT
                variants.vid,
                variants.chrom,
                variants.pos,
                variants.ref,
                variants.alt,
                vals_wide.af,
                vals_wide.dp,
                vals_wide.gt
            FROM
                variants
            INNER JOIN
                vals_wide
            ON 
                variants.vid = vals_wide.vid
            WHERE
                vals_wide.af >= {min_af}
                AND
                vals_wide.dp >= {min_depth}
                AND
                variants.vid NOT IN (
                    SELECT
                        vid
                    FROM
                        filters
                    WHERE
                        filter IN (
                            'map_qual',
                            'slippage',
                            'strand_bias',
                            'weak_evidence',
                            'clustered_events',
                            'base_qual'
                        )
                )
        """)
