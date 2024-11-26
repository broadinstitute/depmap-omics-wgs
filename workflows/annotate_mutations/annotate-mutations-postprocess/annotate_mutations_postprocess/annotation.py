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
                DepMap_ID VARCHAR,
                af FLOAT,
                alt_count UINTEGER,
                am_class VARCHAR,
                am_pathogenicity FLOAT,
                brca1_func_score FLOAT,
                civic_description VARCHAR,
                civic_id VARCHAR,
                civic_score FLOAT,
                dbsnp_rs_id VARCHAR,
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
            CREATE OR REPLACE VIEW transcript_likely_lof_v AS (
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
                    string_agg(transcript_id, ';') AS transcript_likely_lof
                FROM
                    split_to_cols
                WHERE
                    score >= 0.7
                GROUP BY
                    vid
            )
        """)

        db.sql("""
            CREATE OR REPLACE VIEW nmd_v AS (
                WITH exploded AS (
                    SELECT
                        vid,
                        unnest(v_json_arr) AS v_json
                    FROM
                        info
                    WHERE
                        k = 'nmd'
                )
                SELECT
                    vid,
                    list_aggregate(v_json->>'$.*', 'string_agg', '|') AS nmd
                FROM
                    exploded
            )
        """)

        db.sql("""
            CREATE OR REPLACE VIEW hgnc_v AS (
                WITH hgnc_name_exploded AS (
                    SELECT
                        vid,
                        unnest(v_varchar_arr) AS v_varchar
                    FROM
                        info
                    WHERE
                        k = 'hgnc_name'
                ),
                hgnc_name_concat AS (
                    SELECT
                        vid,
                        string_agg(v_varchar, '; ') AS hgnc_name
                    FROM
                        hgnc_name_exploded
                    GROUP BY
                        vid
                ),
                hgnc_group_exploded AS (
                    SELECT
                        vid,
                        unnest(v_varchar_arr) AS v_varchar
                    FROM
                        info
                    WHERE
                        k = 'hgnc_group'
                ),
                hgnc_group_concat AS (
                    SELECT
                        vid,
                        string_agg(v_varchar, '; ') AS hgnc_group
                    FROM
                        hgnc_group_exploded
                    GROUP BY
                        vid
                )
                SELECT
                    coalesce(hgnc_name_concat.vid, hgnc_group_concat.vid) AS vid,
                    hgnc_name_concat.hgnc_name,
                    hgnc_group_concat.hgnc_group
                FROM
                    hgnc_name_concat
                FULL OUTER JOIN
                    hgnc_group_concat
                ON
                    hgnc_name_concat.vid = hgnc_group_concat.vid
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
                                "biotype": "VARCHAR",
                                "clin_sig": "VARCHAR",
                                "consequence": "VARCHAR",
                                "ensp": "VARCHAR",
                                "existing_variation": "VARCHAR",
                                "exon": "VARCHAR",
                                "feature": "VARCHAR",
                                "gene": "VARCHAR",
                                "gnom_ade_af": "FLOAT",
                                "gnom_adg_af": "FLOAT",
                                "hgnc_id": "VARCHAR",
                                "hgvsc": "VARCHAR",
                                "hgvsp": "VARCHAR",
                                "impact": "VARCHAR",
                                "intron": "VARCHAR",
                                "loftool": "FLOAT",
                                "mane_select": "VARCHAR",
                                "pli_gene_value": "FLOAT",
                                "poly_phen": "VARCHAR",
                                "sift": "VARCHAR",
                                "somatic": "VARCHAR",
                                "swissprot": "VARCHAR",
                                "symbol": "VARCHAR",
                                "uniprot_isoform": "VARCHAR",
                                "variant_class": "VARCHAR"
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

        db.sql("DROP VIEW IF EXISTS vep_view")

        db.sql("""
            CREATE OR REPLACE VIEW vals_wide AS (
                SELECT
                    DISTINCT vals.vid,
                    t_ad.ad[1] AS ref_count,
                    t_ad.ad[2] AS alt_count,
                    t_af.af,
                    t_dp.dp,
                    t_gt.gt,
                    t_ps.ps
                FROM
                    vals
                
                LEFT OUTER JOIN (
                    SELECT
                        vid,
                        v_integer_arr AS ad
                    FROM
                        vals
                    WHERE
                        k = 'ad'
                ) t_ad ON vals.vid = t_ad.vid 
                
                LEFT OUTER JOIN (
                    SELECT
                        vid,
                        v_float AS af
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

        db.sql("""
            CREATE OR REPLACE VIEW info_wide AS (
                SELECT
                    DISTINCT info.vid,
                    t_oc_brca1_func_assay_score.oc_brca1_func_assay_score,
                    t_civic_desc.civic_desc,
                    t_civic_id.civic_id,
                    t_civic_score.civic_score,
                    t_rs.rs[1] AS rs,
                    t_gc_prop.gc_prop,
                    list_aggregate(t_mc.mc, 'string_agg', ',') AS mc,
                    t_oc_gtex_gtex_gene.oc_gtex_gtex_gene,
                    t_oc_gwas_catalog_disease.oc_gwas_catalog_disease,
                    t_oc_gwas_catalog_pmid.oc_gwas_catalog_pmid,
                    t_hess_driver.hess_driver,
                    t_hess_signature.hess_signature,
                    t_oc_pharmgkb_id.oc_pharmgkb_id,
                    t_oc_provean_prediction.oc_provean_prediction
                FROM
                    info

                LEFT OUTER JOIN (
                    SELECT
                        vid,
                        v_float AS oc_brca1_func_assay_score
                    FROM
                        info
                    WHERE
                        k = 'oc_brca1_func_assay_score'
                ) t_oc_brca1_func_assay_score ON info.vid = t_oc_brca1_func_assay_score.vid
                
                LEFT OUTER JOIN (
                    SELECT
                        vid,
                        v_varchar AS civic_desc
                    FROM
                        info
                    WHERE
                        k = 'civic_desc'
                ) t_civic_desc ON info.vid = t_civic_desc.vid
                
                LEFT OUTER JOIN (
                    SELECT
                        vid,
                        v_varchar AS civic_id
                    FROM
                        info
                    WHERE
                        k = 'civic_id'
                ) t_civic_id ON info.vid = t_civic_id.vid
                
                LEFT OUTER JOIN (
                    SELECT
                        vid,
                        v_float AS civic_score
                    FROM
                        info
                    WHERE
                        k = 'civic_score'
                ) t_civic_score ON info.vid = t_civic_score.vid
                
                LEFT OUTER JOIN (
                    SELECT
                        vid,
                        v_varchar_arr AS rs
                    FROM
                        info
                    WHERE
                        k = 'rs'
                ) t_rs ON info.vid = t_rs.vid
                
                LEFT OUTER JOIN (
                    SELECT
                        vid,
                        v_float AS gc_prop
                    FROM
                        info
                    WHERE
                        k = 'gc_prop'
                ) t_gc_prop ON info.vid = t_gc_prop.vid
                
                LEFT OUTER JOIN (
                    SELECT
                        vid,
                        v_varchar_arr AS mc
                    FROM
                        info
                    WHERE
                        k = 'mc'
                ) t_mc ON info.vid = t_mc.vid
                
                LEFT OUTER JOIN (
                    SELECT
                        vid,
                        v_varchar AS oc_gtex_gtex_gene
                    FROM
                        info
                    WHERE
                        k = 'oc_gtex_gtex_gene'
                ) t_oc_gtex_gtex_gene ON info.vid = t_oc_gtex_gtex_gene.vid
                
                LEFT OUTER JOIN (
                    SELECT
                        vid,
                        v_varchar AS oc_gwas_catalog_disease
                    FROM
                        info
                    WHERE
                        k = 'oc_gwas_catalog_disease'
                ) t_oc_gwas_catalog_disease ON info.vid = t_oc_gwas_catalog_disease.vid
                
                LEFT OUTER JOIN (
                    SELECT
                        vid,
                        v_varchar AS oc_gwas_catalog_pmid
                    FROM
                        info
                    WHERE
                        k = 'oc_gwas_catalog_pmid'
                ) t_oc_gwas_catalog_pmid ON info.vid = t_oc_gwas_catalog_pmid.vid
                
                LEFT OUTER JOIN (
                    SELECT
                        vid,
                        TRUE AS hess_driver
                    FROM
                        info
                    WHERE
                        k = 'hess'
                ) t_hess_driver ON info.vid = t_hess_driver.vid
                
                LEFT OUTER JOIN (
                    SELECT
                        vid,
                        v_varchar AS hess_signature
                    FROM
                        info
                    WHERE
                        k = 'hess'
                ) t_hess_signature ON info.vid = t_hess_signature.vid
                
                LEFT OUTER JOIN (
                    SELECT
                        vid,
                        v_varchar AS oc_pharmgkb_id
                    FROM
                        info
                    WHERE
                        k = 'oc_pharmgkb_id'
                ) t_oc_pharmgkb_id ON info.vid = t_oc_pharmgkb_id.vid
                
                LEFT OUTER JOIN (
                    SELECT
                        vid,
                        v_varchar AS oc_provean_prediction
                    FROM
                        info
                    WHERE
                        k = 'oc_provean_prediction'
                ) t_oc_provean_prediction ON info.vid = t_oc_provean_prediction.vid
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
                vals_wide.ref_count AS ref_count,
                vals_wide.alt_count AS alt_count,
                vals_wide.af AS af,
                vals_wide.dp AS dp,
                vals_wide.gt AS gt,
                vals_wide.ps AS ps,
                info_wide.oc_brca1_func_assay_score AS brca1_func_score,
                info_wide.civic_desc AS civic_description,
                info_wide.civic_id AS civic_id,
                info_wide.civic_score AS civic_score,
                info_wide.rs AS dbsnp_rs_id,
                info_wide.gc_prop AS gc_content,
                info_wide.mc AS molecular_consequence,
                info_wide.oc_gtex_gtex_gene AS gtex_gene,
                info_wide.oc_gwas_catalog_disease AS gwas_disease,
                info_wide.oc_gwas_catalog_pmid AS gwas_pmid,
                info_wide.hess_driver AS hess_driver,
                info_wide.hess_signature AS hess_signature,
                info_wide.oc_pharmgkb_id AS pharmgkb_id,
                info_wide.oc_provean_prediction AS provean_prediction,
                vep.am_class AS am_class,
                vep.am_pathogenicity AS am_pathogenicity,
                vep.hgvsc AS dna_change,
                vep.feature AS ensembl_feature_id,
                vep.gene AS ensembl_gene_id,
                vep.exon AS exon,
                vep.gnom_ade_af AS gnomade_af,
                vep.gnom_adg_af AS gnomadg_af,
                vep.symbol AS hugo_symbol,
                vep.intron AS intron,
                vep.poly_phen AS polyphen,
                vep.hgvsp AS protein_change,
                vep.sift AS sift,
                vep.uniprot_isoform AS uniprot_id,
                vep.consequence AS variant_info,
                vep.variant_class AS variant_type,
                vep.biotype AS vep_biotype,
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
                transcript_likely_lof_v.transcript_likely_lof AS transcript_likely_lof,
                nmd_v.nmd AS nmd,
                hgnc_v.hgnc_name AS hgnc_name,
                hgnc_v.hgnc_group AS hgnc_family
            FROM
                variants
            INNER JOIN
                vals_wide
            ON 
                variants.vid = vals_wide.vid
            LEFT JOIN
                info_wide
            ON 
                variants.vid = info_wide.vid
            LEFT JOIN
                vep
            ON 
                variants.vid = vep.vid
            LEFT JOIN
                transcript_likely_lof_v
            ON 
                variants.vid = transcript_likely_lof_v.vid
            LEFT JOIN
                nmd_v
            ON 
                variants.vid = nmd_v.vid
            LEFT JOIN
                hgnc_v
            ON 
                variants.vid = hgnc_v.vid
            WHERE
                variants.vid IN (SELECT vid from filtered_vids)
                and
                chrom is not null
        """)
