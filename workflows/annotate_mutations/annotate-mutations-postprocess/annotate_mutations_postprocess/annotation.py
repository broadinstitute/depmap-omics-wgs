import logging
import os
from pathlib import Path

import duckdb


def annotate_vcf(
    db_path: Path,
    parquet_dir_path: Path,
    min_af: float = 0.15,
    min_depth: int = 2,
    max_pop_af: float = 1e-05,
) -> None:
    logging.info("Annotating VCF")

    try:
        os.remove(db_path)
    except OSError:
        pass

    with duckdb.connect(db_path) as db:
        logging.info(f"Reading schema and Parquet files from {parquet_dir_path}")
        db.sql(f"IMPORT DATABASE '{parquet_dir_path}'")

        db.sql("""
            CREATE TABLE IF NOT EXISTS maf (
                chrom VARCHAR,
                pos VARCHAR,
                ref VARCHAR,
                alt VARCHAR,
                af FLOAT,
                alt_count UINTEGER,
                am_class VARCHAR,
                am_pathogenicity FLOAT,
                brca1_func_score FLOAT,
                civic_description VARCHAR,
                civic_id VARCHAR,
                civic_score FLOAT,
                dbsnp_rs_id VARCHAR,
                DepMap_ID VARCHAR,
                dida_id VARCHAR,
                dida_name VARCHAR,
                dna_change VARCHAR,
                dp VARCHAR,
                ensembl_feature_id VARCHAR,
                ensembl_gene_id VARCHAR,
                exon VARCHAR,
                gc_content FLOAT,
                gnomade_af FLOAT,
                gnomadg_af FLOAT,
                gt VARCHAR,
                gtex_gene VARCHAR,
                gwas_disease VARCHAR,
                gwas_pmid VARCHAR,
                hess_driver BOOLEAN,
                hess_signature VARCHAR,
                hgnc_family VARCHAR,
                hgnc_name VARCHAR,
                hugo_symbol VARCHAR,
                intron VARCHAR,
                likely_lof BOOLEAN,
                lof_gene_id VARCHAR,
                lof_gene_name VARCHAR,
                lof_number_of_transcripts_in_gene UINTEGER,
                lof_percent_of_transcripts_affected FLOAT,
                molecular_consequence VARCHAR,
                nmd VARCHAR,
                oncogene_high_impact BOOLEAN,
                pharmgkb_id VARCHAR,
                polyphen VARCHAR,
                protein_change VARCHAR,
                provean_prediction VARCHAR,
                ps VARCHAR,
                ref_count UINTEGER,
                rescue BOOLEAN,
                revel_score FLOAT,
                sift VARCHAR,
                transcript_likely_lof VARCHAR,
                tumor_suppressor_high_impact BOOLEAN,
                uniprot_id VARCHAR,
                variant_info VARCHAR,
                variant_type VARCHAR,
                vep_biotype VARCHAR,
                vep_clin_sig VARCHAR,
                vep_ensp VARCHAR,
                vep_existing_variation VARCHAR,
                vep_hgnc_id VARCHAR,
                vep_impact VARCHAR,
                vep_loftool FLOAT,
                vep_mane_select VARCHAR,
                vep_pli_gene_value FLOAT,
                vep_somatic VARCHAR,
                vep_swissprot VARCHAR
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
            CREATE OR REPLACE VIEW vep_view AS (
                WITH vep_exploded AS (
                    SELECT
                        vid,
                        json_transform(
                            UNNEST(v_json_arr),
                            '{
                                "am_class": "VARCHAR",
                                "am_pathogenicity": "FLOAT",
                                "hgvsc": "VARCHAR",
                                "feature": "VARCHAR",
                                "gene": "VARCHAR",
                                "exon": "VARCHAR",
                                "gnom_ade_af": "FLOAT",
                                "gnom_adg_af": "FLOAT",
                                "symbol": "VARCHAR",
                                "intron": "VARCHAR",
                                "poly_phen": "VARCHAR",
                                "hgvsp": "VARCHAR",
                                "sift": "VARCHAR",
                                "uniprot_isoform": "VARCHAR",
                                "consequence": "VARCHAR",
                                "variant_class": "VARCHAR",
                                "biotype": "VARCHAR",
                                "clin_sig": "VARCHAR",
                                "ensp": "VARCHAR",
                                "existing_variation": "VARCHAR",
                                "hgnc_id": "VARCHAR",
                                "impact": "VARCHAR",
                                "loftool": "FLOAT",
                                "mane_select": "VARCHAR",
                                "pli_gene_value": "FLOAT",
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

        # this is a physical table due to the slowness of the `vep_exploded` CTE
        db.sql("""
            CREATE OR REPLACE TABLE vep
            AS
            SELECT * FROM vep_view
        """)

        db.sql("""
            CREATE OR REPLACE VIEW vals_wide AS (
                SELECT
                    DISTINCT vals.vid,
                    t_ad.ad,
                    t_af.af,
                    t_dp.dp,
                    t_gt.gt,
                    t_ps.ps
                FROM
                    vals
                
                LEFT OUTER JOIN (
                    SELECT
                        vid,
                        v_integer_arr AS ad,
                    FROM
                        vals
                    WHERE
                        k = 'ad'
                ) t_ad ON vals.vid = t_ad.vid 
                
                LEFT OUTER JOIN (
                    SELECT
                        vid,
                        v_float AS af,
                    FROM
                        vals
                    WHERE
                        k = 'af'
                ) t_af ON vals.vid = t_af.vid 
                
                LEFT OUTER JOIN (
                    SELECT
                        vid,
                        v_integer AS dp
                    FROM
                        vals
                    WHERE
                        k = 'dp'
                ) t_dp ON vals.vid = t_dp.vid
                
                LEFT OUTER JOIN (
                    SELECT
                        vid,
                        v_varchar AS gt
                    FROM
                        vals
                    WHERE
                        k = 'gt'
                ) t_gt ON vals.vid = t_gt.vid 
                
                LEFT OUTER JOIN (
                    SELECT
                        vid,
                        v_integer AS ps
                    FROM
                        vals
                    WHERE
                        k = 'ps'
                ) t_ps ON vals.vid = t_ps.vid 
            )
        """)

        info_col_name_map = {
            "civic_desc": "civic_description",
            "civic_score": "civic_score",
            "oc_revel_score": "revel_score",
            "oc_pharmgkb_id": "pharmgkb_id",
            "gc_prop": "gc_prop",
        }

        maf_info_cols = db.execute(
            """
            SELECT
                id_snake,
                v_col_name
            FROM
                val_info_types
            WHERE
                kind = 'info'
                AND
                id_snake IN (
                    SELECT DISTINCT k FROM info
                )
                AND
                ('rs' IN $maf_info_cols)
        """,
            {"maf_info_cols": ["dp", "rs"]},
        ).df()

        relevent_info = [
            "civic_desc",
            "civic_score",
            "oc_revel_score",
            "oc_pharmgkb_id",
            "gc_prop",
        ]

        relevent_info_clause = ", ".join(f"'{x}'" for x in relevent_info)

        db.execute(
            """
            select * from info where k IN $relevent_info
        """,
            {"relevent_info": relevent_info},
        )

        pivot_ctes = [
            f"""
            info_{c} AS (
                PIVOT (
                    SELECT
                        vid,
                        k,
                        {c}
                    FROM
                        info
                    WHERE
                        --k IN ({relevent_info})
                        --AND
                        {c} IS NOT NULL
                ) ON k USING first({c})
            )
        """
            for c in info_cols
        ]

        join_clauses = [
            f"""
            FULL OUTER JOIN
                info_{c}
            ON
                info_{info_cols[0]}.vid = info_{c}.vid
        """
            for c in info_cols[1:]
        ]

        # this is a physical table because you can't PIVOT this way inside views
        db.sql(f"""
            WITH {', '.join(pivot_ctes)}
            SELECT * FROM info_{info_cols[0]}
            {'\n'.join(join_clauses)}
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
            CREATE OR REPLACE VIEW filtered_vids AS (
                SELECT
                    vid
                FROM
                    variants
                WHERE
                    -- quality checks
                    vid IN (
                        SELECT
                            vid
                        FROM
                            vals_wide
                        WHERE
                            vals_wide.af >= {min_af}
                            AND
                            vals_wide.dp >= {min_depth}
                    )
                    AND
                    vid NOT IN (
                        SELECT
                            vid
                        FROM
                            filters
                        WHERE
                            filter IN (
                                'multiallelic',
                                'map_qual',
                                'slippage',
                                'strand_bias',
                                'weak_evidence',
                                'clustered_events',
                                'base_qual'
                            )
                    )
                    AND
                    vid IN (
                        SELECT
                            vid
                        FROM
                            vep
                        WHERE
                            -- important splice event
                            (
                                contains(consequence, 'splice')
                                AND
                                impact IN ('HIGH', 'MODERATE')
                            )
                            OR
                            -- protein sequence has changed
                            (
                                hgvsp IS NOT NULL
                                AND
                                NOT ends_with(hgvsp, '=')
                            )
                    )
                    AND
                    -- not in a segmental duplication nor repeatmasker region
                    vid NOT IN (
                        SELECT
                            vid
                        FROM
                            info
                        WHERE
                            k in ('segdup', 'rm', 'pon')
                            AND
                            v_boolean
                    )
                    AND
                    -- max population prevalence per gnomAD
                    vid IN (
                        SELECT
                            vid
                        FROM
                            vep
                        WHERE
                            coalesce(gnom_ade_af, 0) <= {max_pop_af}
                            AND
                            coalesce(gnom_adg_af, 0) <= {max_pop_af}
                    )
            )
        """)

        db.sql("""
            SELECT
                variants.vid,
                variants.chrom AS chrom,
                variants.pos AS pos,
                variants.ref AS ref,
                variants.alt AS alt,
                vals_wide.ad[1] AS ref_count,
                vals_wide.ad[2] AS alt_count,
                vals_wide.af AS af,
                vals_wide.dp AS dp,
                vals_wide.gt AS gt,
                brca1.oc_brca1_func_assay_score AS brca1_func_score,
                vep.clin_sig AS vep_clin_sig,
                vep.ensp AS vep_ensp,
                vep.existing_variation AS vep_existing_variation,
                vep.hgnc_id AS vep_hgnc_id,
                vep.impact AS vep_impact,
                vep.loftool AS vep_loftool,
                vep.mane_select AS vep_mane_select,
                vep.pli_gene_value AS vep_pli_gene_value,
                vep.somatic AS vep_somatic,
                vep.swissprot AS vep_swissprot,
                vep.exon AS exon,
                vep.intron AS intron,
                vep.poly_phen AS polyphen,
                vep.am_class AS am_class,
                vep.am_pathogenicity AS am_pathogenicity,
                list_aggregate(info_wide.mc, 'string_agg', ',') AS molecular_consequence,
                NULL AS civic_id,
                info_wide.civic_desc AS civic_description,
                info_wide.civic_score AS civic_score,
                info_wide.oc_revel_score AS revel_score,
                info_wide.oc_pharmgkb_id AS pharmgkb_id,
                info_wide.gc_prop AS gc_prop,
            FROM
                variants
            INNER JOIN
                vals_wide
            ON 
                variants.vid = vals_wide.vid
            LEFT JOIN
                brca1
            ON 
                variants.vid = brca1.vid
            LEFT JOIN
                vep
            ON 
                variants.vid = vep.vid
            LEFT JOIN
                revel
            ON 
                variants.vid = revel.vid
            LEFT JOIN
                info_wide
            ON 
                variants.vid = info_wide.vid
            WHERE
                variants.vid IN (SELECT vid from filtered_vids)
                and
                chrom is not null
        """)
