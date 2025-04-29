import logging
import os
from pathlib import Path

import duckdb
import pandas as pd


def convert_duckdb_to_maf(
    db_path: Path,
    parquet_dir_path: Path,
    out_file_path: Path,
    min_af: float = 0.15,
    min_depth: int = 5,
    max_pop_af: float = 1e-05,
    max_brca1_func_assay_score: float = -1.328,
) -> None:
    """
    Convert a vcf-to-duckdb DuckDB database containing a vcf-to-duckdb to a MAF-like
    wide data frame.

    Reads schema and Parquet files from the specified directory, creates necessary views
    and tables in the database, filters germline and low-quality variants, and exports
    the resulting somatic/rescued variants to a Parquet file.

    Arguments:
        db_path: Path to the DuckDB database file
        parquet_dir_path: Directory containing Parquet files to import
        out_file_path: Output file path for the MAF file (in Parquet format)
        min_af: Minimum allele frequency threshold for variant filtering
        min_depth: Minimum read depth threshold for variant filtering
        max_pop_af: Maximum population allele frequency for variant filtering
        max_brca1_func_assay_score: Maximum BRCA1 functional assay score for rescuing
    """

    logging.info("Annotating VCF")

    try:
        os.remove(db_path)
    except OSError:
        pass

    with duckdb.connect(db_path) as db:
        logging.info(f"Reading schema and Parquet files from {parquet_dir_path}")
        db.sql(f"IMPORT DATABASE '{parquet_dir_path}'")

        make_views(db, min_af, min_depth, max_pop_af, max_brca1_func_assay_score)

        somatic_variants = get_somatic_variants_as_df(db)
        somatic_variants.to_parquet(out_file_path, index=False)


def make_views(
    db: duckdb.DuckDBPyConnection,
    min_af: float = 0.15,
    min_depth: int = 5,
    max_pop_af: float = 1e-05,
    max_brca1_func_assay_score: float = -1.328,
):
    """
    Create all necessary views and tables in the DuckDB database.

    Sets up the database structure by creating views for filtering variants based on
    quality criteria, extracting information from the vals_info table, and preparing
    data for export to MAF format.

    Arguments:
        db: DuckDB database connection
        min_af: Minimum allele frequency threshold for variant filtering
        min_depth: Minimum read depth threshold for variant filtering
        max_pop_af: Maximum population allele frequency for variant filtering
        max_brca1_func_assay_score: Maximum BRCA1 functional assay score for rescuing
    """

    # separate vals_info table into two views
    make_vals_info_views(db)

    # convert vals to wide view
    make_vals_wide_view(db)

    # make views for all filters and global filter for high-quality variants
    make_filters_view(db)
    make_quality_vids_view(db, min_af, min_depth)

    # all views/queries after this point should filter on quality_vids
    make_info_wide_view(db)
    make_transcript_likely_lof_view(db)
    make_nmd_view(db)
    make_lof_view(db)
    make_hgnc_view(db)
    make_vep_table(db)
    make_oncogene_tsg_view(db)
    make_rescues_view(db, max_brca1_func_assay_score)
    make_filtered_vids_view(db, max_pop_af)
    make_somatic_variants_table(db)


def make_vals_info_views(db: duckdb.DuckDBPyConnection) -> None:
    """
    Create views that separate the vals_info table into 'vals' and 'info' views.

    Splits the vals_info table into two separate views based on the 'kind' column,
    creating a 'vals' view for variant values and an 'info' view for variant information.

    Arguments:
        db: DuckDB database connection
    """

    logging.info("Making vals and info views")

    db.sql("""
        CREATE OR REPLACE VIEW vals AS (
            SELECT
                *
            FROM
                vals_info
            WHERE
                kind = 'val'
        )
    """)

    db.sql("""
        CREATE OR REPLACE VIEW info AS (
            SELECT
                *
            FROM
                vals_info
            WHERE
                kind = 'info'
        )
    """)


def make_vals_wide_view(db: duckdb.DuckDBPyConnection) -> None:
    """
    Create a wide view of variant values from the 'vals' view.

    Transforms the 'vals' view from a long format to a wide format, with columns for
    reference count, alternate count, allele frequency, depth, genotype, and phase set.

    Arguments:
        db: DuckDB database connection
    """

    logging.info("Making vals_wide view")

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


def make_filters_view(db: duckdb.DuckDBPyConnection) -> None:
    """
    Create a view of variant filters from the 'variants' and 'info' tables.

    Combines filters from the 'variants' table and allele-specific filters from the
    'info' table into a single view, with one row per variant-filter combination.

    Arguments:
        db: DuckDB database connection
    """

    logging.info("Making filters view")

    db.sql("""
        CREATE OR REPLACE VIEW filters AS (
            WITH variant_filters AS (
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
                variant_filters
            UNION
            SELECT
                vid,
                lower(filter) AS filter
            FROM
                as_filters
        )
    """)


def make_quality_vids_view(
    db: duckdb.DuckDBPyConnection, min_af: float, min_depth: int
) -> None:
    """
    Create a view of high-quality variant IDs based on filtering criteria.

    Filters variants based on allele frequency, depth, and absence of certain filters
    like multiallelic, map_qual, slippage, strand_bias, weak_evidence, and base_qual.

    Arguments:
        db: DuckDB database connection
        min_af: Minimum allele frequency threshold
        min_depth: Minimum read depth threshold
    """

    logging.info("Making quality_vids view")

    db.sql(f"""
        CREATE OR REPLACE VIEW quality_vids AS (
            SELECT
                vid
            FROM
                variants
            WHERE
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
                            'base_qual'
                        )
                )
        )
    """)


def make_info_wide_view(db: duckdb.DuckDBPyConnection) -> None:
    """
    Create a wide view of variant information from the 'info' table.

    Transforms the 'info' table from a long format to a wide format, with columns for
    various annotations like BRCA1 functional assay score, CIViC information, COSMIC tier,
    dbSNP ID, and other variant annotations.

    Arguments:
        db: DuckDB database connection
    """

    logging.info("Making info_wide view")

    db.sql("""
        CREATE OR REPLACE VIEW info_wide AS (
            SELECT
                DISTINCT info.vid,
                t_oc_brca1_func_assay_score.oc_brca1_func_assay_score,
                t_civic_desc.civic_desc,
                t_civic_id.civic_id,
                t_civic_score.civic_score,
                t_cosmic_tier.cosmic_tier,
                t_rs.rs[1] AS rs,
                t_gc_prop.gc_prop,
                list_aggregate(t_mc.mc, 'string_agg', ',') AS mc,
                t_oc_gtex_gtex_gene.oc_gtex_gtex_gene,
                t_oc_gwas_catalog_disease.oc_gwas_catalog_disease,
                t_oc_gwas_catalog_pmid.oc_gwas_catalog_pmid,
                t_hess.hess_driver,
                t_hess.hess_signature,
                t_oc_pharmgkb_id.oc_pharmgkb_id,
                t_oc_provean_prediction.oc_provean_prediction,
                t_oc_revel_score.oc_revel_score,
                t_oncokb_hotspot.oncokb_hotspot,
                t_oncokb_muteff.oncokb_muteff,
                t_oncokb_oncogenic.oncokb_oncogenic
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
            ) t_oc_brca1_func_assay_score ON
                info.vid = t_oc_brca1_func_assay_score.vid

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
                    v_integer AS civic_id
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
                    v_integer AS cosmic_tier
                FROM
                    info
                WHERE
                    k = 'cmc_tier'
            ) t_cosmic_tier ON info.vid = t_cosmic_tier.vid

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
                    v_varchar AS hess_signature,
                    TRUE AS hess_driver
                FROM
                    info
                WHERE
                    k = 'hess'
            ) t_hess ON info.vid = t_hess.vid

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

            LEFT OUTER JOIN (
                SELECT
                    vid,
                    v_float AS oc_revel_score
                FROM
                    info
                WHERE
                    k = 'oc_revel_score'
            ) t_oc_revel_score ON info.vid = t_oc_revel_score.vid

            LEFT OUTER JOIN (
                SELECT
                    vid,
                    v_varchar AS oncokb_muteff
                FROM
                    info
                WHERE
                    k = 'oncokb_muteff'
            ) t_oncokb_muteff ON info.vid = t_oncokb_muteff.vid

            LEFT OUTER JOIN (
                SELECT
                    vid,
                    v_boolean AS oncokb_hotspot
                FROM
                    info
                WHERE
                    k = 'oncokb_hotspot'
            ) t_oncokb_hotspot ON info.vid = t_oncokb_hotspot.vid

            LEFT OUTER JOIN (
                SELECT
                    vid,
                    v_varchar AS oncokb_oncogenic
                FROM
                    info
                WHERE
                    k = 'oncokb_oncogenic'
            ) t_oncokb_oncogenic ON info.vid = t_oncokb_oncogenic.vid

            WHERE info.vid IN (SELECT vid FROM quality_vids)
        )
    """)


def make_transcript_likely_lof_view(db: duckdb.DuckDBPyConnection) -> None:
    """
    Create a view of transcripts likely to have loss-of-function (LoF) variants.

    Processes the 'oc_revel_all' field from the 'info' table to identify transcripts
    with a REVEL score >= 0.7, which indicates likely loss-of-function.

    Arguments:
        db: DuckDB database connection
    """

    logging.info("Making transcript_likely_lof_v view")

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
                    AND
                    vid IN (SELECT vid FROM quality_vids)
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


def make_nmd_view(db: duckdb.DuckDBPyConnection) -> None:
    """
    Create a view of nonsense-mediated decay (NMD) information.

    Processes the 'nmd' field from the 'info' table to extract information about
    nonsense-mediated decay for each variant.

    Arguments:
        db: DuckDB database connection
    """

    logging.info("Making nmd_v view")

    db.sql("""
        CREATE OR REPLACE VIEW nmd_v AS (
            WITH exploded AS (
                SELECT
                    vid,
                    v_json_arr[1] AS v_json -- SnpEff actually puts only one value here
                FROM
                    info
                WHERE
                    k = 'nmd'
                    AND
                    vid IN (SELECT vid FROM quality_vids)
            )
            SELECT
                vid,
                list_aggregate(v_json->>'$.*', 'string_agg', '|') AS nmd
            FROM
                exploded
        )
    """)


def make_lof_view(db: duckdb.DuckDBPyConnection) -> None:
    """
    Create a view of loss-of-function (LoF) information.

    Processes the 'lof' field from the 'info' table to extract information about
    loss-of-function effects, including gene name, gene ID, number of transcripts,
    and proportion of transcripts affected.

    Arguments:
        db: DuckDB database connection
    """

    logging.info("Making lof_v view")

    db.sql("""
        CREATE OR REPLACE VIEW lof_v AS (
            WITH exploded AS (
                SELECT
                    vid,
                    v_json_arr[1] AS v_json -- SnpEff actually puts only one value here
                FROM
                    info
                WHERE
                    k = 'lof'
                    AND
                    vid IN (SELECT vid FROM quality_vids)
            )
            SELECT
                vid,
                v_json->>'$.gene_name' AS gene_name,
                v_json->>'$.gene_id' AS gene_id,
                (v_json->>'$.number_of_transcripts_in_gene')::USMALLINT
                    AS number_of_transcripts_in_gene,
                (v_json->>'$.percent_of_transcripts_affected')::FLOAT
                    AS prop_of_transcripts_affected
            FROM
                exploded
        )
    """)


def make_hgnc_view(db: duckdb.DuckDBPyConnection) -> None:
    """
    Create a view of HGNC (HUGO Gene Nomenclature Committee) information.

    Processes the 'hgnc_name' and 'hgnc_group' fields from the 'info' table to extract
    information about gene names and gene families.

    Arguments:
        db: DuckDB database connection
    """

    logging.info("Making hgnc_v view")

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
                    AND
                    vid IN (SELECT vid FROM quality_vids)
            ),
            hgnc_name_concat AS (
                SELECT
                    vid,
                    string_agg(v_varchar, ';') AS hgnc_name
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
                    AND
                    vid IN (SELECT vid FROM quality_vids)
            ),
            hgnc_group_concat AS (
                SELECT
                    vid,
                    string_agg(v_varchar, ';') AS hgnc_group
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


def make_vep_table(db: duckdb.DuckDBPyConnection) -> None:
    """
    Create a table of Variant Effect Predictor (VEP) annotations.

    Processes the 'csq' field from the 'info' table to extract VEP annotations for each
    variant. Creates a physical table rather than a view for performance reasons.

    Arguments:
        db: DuckDB database connection
    """

    logging.info("Making vep table")

    # this is a physical table due to the relative slowness of the `vep_exploded` CTE
    # and the fact it gets queried several times later
    db.sql("""
        CREATE TABLE IF NOT EXISTS vep
        AS
        SELECT * FROM (
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
                    AND
                    vid IN (SELECT vid FROM quality_vids)
            )
            SELECT
                vid,
                csq.am_class,
                csq.am_pathogenicity,
                csq.biotype,
                csq.clin_sig,
                csq.consequence,
                csq.ensp,
                csq.existing_variation,
                csq.exon,
                csq.feature,
                csq.gene,
                csq.gnom_ade_af,
                csq.gnom_adg_af,
                csq.hgnc_id,
                csq.hgvsc,
                url_decode(csq.hgvsp) AS hgvsp,
                csq.impact,
                csq.intron,
                csq.loftool,
                csq.mane_select,
                csq.pli_gene_value,
                csq.poly_phen,
                csq.sift,
                csq.somatic,
                csq.swissprot,
                csq.symbol,
                csq.uniprot_isoform,
                csq.variant_class
            FROM
                vep_exploded
        )
    """)


def make_oncogene_tsg_view(db: duckdb.DuckDBPyConnection) -> None:
    """
    Create a view of oncogene and tumor suppressor gene (TSG) information.

    Identifies variants in oncogenes and tumor suppressor genes with high-impact effects,
    based on the 'oncogene' and 'tsg' fields from the 'info' table and the 'impact'
    field from the 'vep' table.

    Arguments:
        db: DuckDB database connection
    """

    logging.info("Making oncogene_tsg view")

    db.sql("""
        CREATE OR REPLACE VIEW oncogene_tsg AS (
            WITH oncogene_v AS (
                SELECT
                    vid,
                    TRUE AS oncogene_high_impact
                FROM
                    info
                WHERE
                    k = 'oncogene'
                    AND
                    vid IN (SELECT vid FROM vep WHERE impact = 'HIGH')
            ),
            tsg_v AS (
                SELECT
                    vid,
                    TRUE AS tumor_suppressor_high_impact
                FROM
                    info
                WHERE
                    k = 'tsg'
                    AND
                    vid IN (SELECT vid FROM vep WHERE impact = 'HIGH')
            )
            SELECT
                coalesce(oncogene_v.vid, tsg_v.vid) AS vid,
                oncogene_v.oncogene_high_impact,
                tsg_v.tumor_suppressor_high_impact
            FROM
                oncogene_v
            FULL OUTER JOIN
                tsg_v
            ON
                oncogene_v.vid = tsg_v.vid
        )
    """)


def make_rescues_view(
    db: duckdb.DuckDBPyConnection, max_brca1_func_assay_score: float
) -> None:
    """
    Create a view of rescued variants based on various criteria.

    Identifies variants that should be rescued despite not meeting standard quality
    criteria, based on factors like oncogenic effect, hotspot status, COSMIC tier,
    BRCA1 functional assay score, oncogene/TSG status, and specific gene-position
    combinations (e.g., TERT promoter, MET exon 14 skipping).

    Arguments:
        db: DuckDB database connection
        max_brca1_func_assay_score: Maximum BRCA1 functional assay score for rescuing
    """

    logging.info("Making rescues view")

    db.sql(f"""
        CREATE OR REPLACE VIEW rescues AS (
            SELECT
                DISTINCT vid,
                TRUE AS rescued
            FROM
                info
            WHERE
                vid IN (SELECT vid FROM quality_vids)
                AND (
                    (
                        k = 'oncokb_muteff'
                        AND
                        v_varchar IN ('Loss-of-function', 'Gain-of-function')
                    )
                    OR
                    (
                        k = 'oncokb_oncogenic'
                        AND
                        v_varchar = 'Oncogenic'
                    )
                    OR
                    (
                        k = 'oncokb_hotspot'
                        AND
                        v_boolean
                    )
                    OR
                    (
                        k = 'cmc_tier'
                        AND
                        v_integer = 1
                    )
                    OR
                    (
                        k = 'oc_brca1_func_assay_score'
                        AND
                        v_float <= {max_brca1_func_assay_score}
                    )
                    OR
                    vid IN (
                        SELECT
                            vid
                        FROM
                            oncogene_tsg
                        WHERE
                            oncogene_high_impact OR tumor_suppressor_high_impact
                    )
                    OR
                    vid IN (
                        SELECT
                            vid
                        FROM
                            info
                        WHERE
                            k = 'hess'
                    )
                    OR
                    vid IN (
                        SELECT
                            variants.vid
                        FROM
                            variants
                        INNER JOIN
                            vep
                        ON
                            variants.vid = vep.vid
                        WHERE
                            variants.chrom = 'chr5'
                            AND
                            variants.pos BETWEEN 1295054 AND 1295365
                            AND
                            vep.symbol = 'TERT'
                    )
                    OR
                    vid IN (
                        SELECT
                            variants.vid
                        FROM
                            variants
                        INNER JOIN
                            vep
                        ON
                            variants.vid = vep.vid
                        WHERE
                            variants.chrom = 'chr7'
                            AND
                            variants.pos BETWEEN 116771825 AND 116771840
                            AND
                            vep.symbol = 'MET'
                    )
                )
            )
    """)


def make_filtered_vids_view(db: duckdb.DuckDBPyConnection, max_pop_af: float) -> None:
    """
    Create a view of filtered variant IDs based on various criteria.

    Identifies variants that pass quality filters and either are rescued or meet criteria
    for being valid somatic alterations (splice events, protein changes) while not being
    in clustered events, segmental duplications, or repeat regions, and having low
    population allele frequency.

    Arguments:
        db: DuckDB database connection
        max_pop_af: Maximum population allele frequency
    """

    logging.info("Making filtered_vids view")

    db.sql(f"""
        CREATE OR REPLACE VIEW filtered_vids AS (
            SELECT
                vid
            FROM
                variants
            WHERE
                vid IN (SELECT vid FROM quality_vids)
                AND (
                    -- vid is rescued
                    vid IN (SELECT vid FROM rescues)
                    -- vid is a valid somatic alteration
                    OR (
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
                        -- not part of a clustered event according to mutect2
                        vid NOT IN (
                            SELECT
                                vid
                            FROM
                                filters
                            WHERE
                                filter = 'clustered_event'
                        )
                        AND
                        -- not in a segmental duplication nor repeatmasker region
                        vid NOT IN (
                            SELECT
                                vid
                            FROM
                                info
                            WHERE
                                k IN ('segdup', 'rm')
                                AND
                                v_boolean
                        )
                        AND
                        -- below max population prevalence per gnomAD
                        vid IN (
                            SELECT
                                vid
                            FROM
                                vep
                            WHERE
                                coalesce(gnom_ade_af, 0) <= {max_pop_af}
                                AND
                                coalesce(gnom_adg_af, 0) <= {max_pop_af}
                                AND
                                vid NOT IN (
                                    SELECT
                                        vid
                                    FROM
                                        info
                                    WHERE
                                        k  = 'pon'
                                        AND
                                        v_boolean
                                )
                        )
                    )
                )
        )
    """)


def make_somatic_variants_table(db: duckdb.DuckDBPyConnection) -> None:
    """
    Create a table of somatic variants with all annotations.

    Combines information from all previously created views and tables to create a
    comprehensive table of somatic variants with all annotations. Checks for duplicate
    variants that might indicate a fan trap in the SQL queries.

    Arguments:
        db: DuckDB database connection
    """

    logging.info("Making somatic_variants table")

    db.sql("""
        CREATE TABLE IF NOT EXISTS somatic_variants (
            chrom VARCHAR,
            pos UINTEGER,
            ref VARCHAR,
            alt VARCHAR,
            af DECIMAL(),
            alt_count UINTEGER,
            am_class VARCHAR,
            am_pathogenicity FLOAT,
            brca1_func_score FLOAT,
            civic_description VARCHAR,
            civic_id VARCHAR,
            civic_score FLOAT,
            cosmic_tier UINTEGER,
            dbsnp_rs_id VARCHAR,
            dna_change VARCHAR,
            dp USMALLINT,
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
            lof_gene_id VARCHAR,
            lof_gene_name VARCHAR,
            lof_number_of_transcripts_in_gene UINTEGER,
            lof_prop_of_transcripts_affected FLOAT,
            molecular_consequence VARCHAR,
            nmd VARCHAR,
            oncogene_high_impact BOOLEAN,
            oncokb_effect VARCHAR,
            oncokb_hotspot BOOLEAN,
            oncokb_oncogenic VARCHAR,
            pharmgkb_id VARCHAR,
            polyphen VARCHAR,
            protein_change VARCHAR,
            provean_prediction VARCHAR,
            ps UINTEGER,
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
        INSERT INTO
            somatic_variants
        BY NAME (
            SELECT
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
                info_wide.cosmic_tier AS cosmic_tier,
                info_wide.rs AS dbsnp_rs_id,
                info_wide.gc_prop AS gc_content,
                info_wide.mc AS molecular_consequence,
                info_wide.oc_gtex_gtex_gene AS gtex_gene,
                info_wide.oc_gwas_catalog_disease AS gwas_disease,
                info_wide.oc_gwas_catalog_pmid AS gwas_pmid,
                info_wide.oc_pharmgkb_id AS pharmgkb_id,
                info_wide.oc_provean_prediction AS provean_prediction,
                info_wide.oc_revel_score AS revel_score,
                info_wide.oncokb_muteff AS oncokb_effect,
                coalesce(info_wide.oncokb_hotspot, FALSE) AS oncokb_hotspot,
                info_wide.oncokb_oncogenic AS oncokb_oncogenic,
                coalesce(info_wide.hess_driver, FALSE) AS hess_driver,
                info_wide.hess_signature AS hess_signature,
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
                hgnc_v.hgnc_group AS hgnc_family,
                lof_v.gene_id AS lof_gene_id,
                lof_v.gene_name AS lof_gene_name,
                lof_v.number_of_transcripts_in_gene AS
                    lof_number_of_transcripts_in_gene,
                lof_v.prop_of_transcripts_affected AS
                    lof_prop_of_transcripts_affected,
                coalesce(oncogene_tsg.oncogene_high_impact, FALSE) AS
                    oncogene_high_impact,
                coalesce(oncogene_tsg.tumor_suppressor_high_impact, FALSE) AS
                    tumor_suppressor_high_impact,
                coalesce(rescues.rescued, FALSE) AS rescue
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
            LEFT JOIN
                lof_v
            ON
                variants.vid = lof_v.vid
            LEFT JOIN
                oncogene_tsg
            ON
                variants.vid = oncogene_tsg.vid
            LEFT JOIN
                rescues
            ON
                variants.vid = rescues.vid
            WHERE
                variants.vid IN (SELECT vid from filtered_vids)
            ORDER BY
                replace(replace(chrom[4:], 'X', '23'), 'Y', '24')::INTEGER,
                pos,
                alt
        )
    """)

    # if we end up with duplicate variants, we must have introduced a fan trap due to a
    # missing `group by` somewhere
    dups = db.sql("""
        SELECT
            chrom,
            pos,
            ref,
            alt
        FROM
            somatic_variants
        GROUP BY
            chrom,
            pos,
            ref,
            alt
        HAVING
            count(*) > 1
    """)

    if len(dups) > 0:
        raise ValueError(f"{len(dups)} duplicate variants in final result set:\n{dups}")


def get_somatic_variants_as_df(db: duckdb.DuckDBPyConnection) -> pd.DataFrame:
    """
    Retrieve somatic variants from the database as a pandas DataFrame.

    Queries the 'somatic_variants' table and converts the results to a pandas DataFrame
    with appropriate data types for each column.

    Arguments:
        db: DuckDB database connection

    Returns:
        A pandas DataFrame containing somatic variants with all annotations
    """

    return (
        db.table("somatic_variants")
        .df()
        .astype(
            {
                "chrom": "string",
                "pos": "UInt32",
                "ref": "string",
                "alt": "string",
                "af": "Float32",
                "alt_count": "UInt32",
                "am_class": "string",
                "am_pathogenicity": "Float32",
                "brca1_func_score": "Float32",
                "civic_description": "string",
                "civic_id": "string",
                "civic_score": "Float32",
                "cosmic_tier": "UInt32",
                "dbsnp_rs_id": "string",
                "dna_change": "string",
                "dp": "UInt32",
                "ensembl_feature_id": "string",
                "ensembl_gene_id": "string",
                "exon": "string",
                "gc_content": "Float32",
                "gnomade_af": "Float32",
                "gnomadg_af": "Float32",
                "gt": "string",
                "gtex_gene": "string",
                "gwas_disease": "string",
                "gwas_pmid": "string",
                "hess_driver": "boolean",
                "hess_signature": "string",
                "hgnc_family": "string",
                "hgnc_name": "string",
                "hugo_symbol": "string",
                "intron": "string",
                "lof_gene_id": "string",
                "lof_gene_name": "string",
                "lof_number_of_transcripts_in_gene": "UInt32",
                "lof_prop_of_transcripts_affected": "Float32",
                "molecular_consequence": "string",
                "nmd": "string",
                "oncogene_high_impact": "boolean",
                "oncokb_effect": "string",
                "oncokb_hotspot": "boolean",
                "oncokb_oncogenic": "string",
                "pharmgkb_id": "string",
                "polyphen": "string",
                "protein_change": "string",
                "provean_prediction": "string",
                "ps": "UInt32",
                "ref_count": "UInt32",
                "rescue": "boolean",
                "revel_score": "Float32",
                "sift": "string",
                "transcript_likely_lof": "string",
                "tumor_suppressor_high_impact": "boolean",
                "uniprot_id": "string",
                "variant_info": "string",
                "variant_type": "string",
                "vep_biotype": "string",
                "vep_clin_sig": "string",
                "vep_ensp": "string",
                "vep_existing_variation": "string",
                "vep_hgnc_id": "string",
                "vep_impact": "string",
                "vep_loftool": "Float32",
                "vep_mane_select": "string",
                "vep_pli_gene_value": "Float32",
                "vep_somatic": "string",
                "vep_swissprot": "string",
            }
        )
    )
