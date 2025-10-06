import logging
import os
from math import ceil
from pathlib import Path

import duckdb
import pandas as pd


def do_select_somatic_variants(
    db_path: Path,
    parquet_dir_path: Path,
    somatic_variants_out_file_path: Path,
    variants_enriched_out_file_path: Path | None,
    mut_sig_out_file_path: Path | None,
    db_tmp_dir_path: Path,
    sample_id: str,
    min_af: float = 0.15,
    min_depth: int = 5,
    max_pop_af: float = 1e-05,
    max_brca1_func_assay_score: float = -1.328,
    batch_size: int = 1000000,
) -> None:
    """
    Convert a vcf-to-duckdb DuckDB database to a MAF-like wide table, exported as a
    parquet file.

    Reads schema and Parquet files from the specified directory, creates necessary views
    and tables in the database, filters variants based on quality criteria, and exports
    the resulting somatic/rescued variants to a Parquet file.

    :param db_path: Path to the DuckDB database file
    :param parquet_dir_path: Directory containing Parquet files to import
    :param somatic_variants_out_file_path: Output file path for the somatic variants
    file (in Parquet format)
    :param variants_enriched_out_file_path: Output file path for a wide file of all
    variants (in Parquet format)
    :param mut_sig_out_file_path: Output file path for a narrower version of
    variants_enriched, e.g. for mutational signatures (in Parquet format)
    :param db_tmp_dir_path: path to a temporary directory to use for DuckDB offloading
    :param sample_id: the sample ID
    :param min_af: Minimum allele frequency threshold for variant filtering
    :param min_depth: Minimum read depth threshold for variant filtering
    :param max_pop_af: Maximum population allele frequency for variant filtering
    :param max_brca1_func_assay_score: Maximum BRCA1 functional assay score for variants
    :param batch_size: number of variants to process at a time
    """

    try:
        os.remove(db_path)
    except OSError:
        pass

    with duckdb.connect(db_path) as db:
        # set some misc. DuckDB settings
        db.execute(f"PRAGMA temp_directory='{db_tmp_dir_path.name}'")
        db.execute("SET preserve_insertion_order = false;")

        logging.info(f"Reading schema and Parquet files from {parquet_dir_path}")
        db.sql(f"IMPORT DATABASE '{parquet_dir_path}'")

        # make views and tables
        make_views_and_tables(
            db,
            sample_id,
            min_af,
            min_depth,
            max_pop_af,
            max_brca1_func_assay_score,
            batch_size,
        )

        if variants_enriched_out_file_path is not None:
            # write wide version of all quality variants to parquet
            r = db.table("variants_enriched").count("*").fetchone()
            n_variants = r[0]  # pyright: ignore
            logging.info(
                f"Writing {n_variants}"
                f" enriched variants to {variants_enriched_out_file_path}"
            )
            db.sql(f"""
                COPY
                    variants_enriched
                TO
                    '{variants_enriched_out_file_path}'
                (
                    OVERWRITE,
                    FORMAT parquet,
                    COMPRESSION zstd,
                    COMPRESSION_LEVEL 22
                )
            """)

        if mut_sig_out_file_path is not None:
            # write narrower version of variants for mutational signatures to parquet
            r = db.table("mut_sig").count("*").fetchone()
            n_variants = r[0]  # pyright: ignore
            logging.info(
                f"Writing {n_variants}"
                f" mutational signature variants to {mut_sig_out_file_path}"
            )
            db.sql(f"""
                COPY
                    mut_sig
                TO
                    '{mut_sig_out_file_path}'
                (
                    OVERWRITE,
                    FORMAT parquet,
                    COMPRESSION zstd,
                    COMPRESSION_LEVEL 22
                )
            """)

        # write the somatic variants to parquet (via data frame to ensure
        # backwards-compatible dtypes)
        somatic_variants = get_somatic_variants_as_df(db)
        logging.info(
            f"Writing {somatic_variants.shape[0]} somatic variants"
            f" to {somatic_variants_out_file_path}"
        )
        somatic_variants.to_parquet(somatic_variants_out_file_path, index=False)


def make_views_and_tables(
    db: duckdb.DuckDBPyConnection,
    sample_id: str,
    min_af: float = 0.15,
    min_depth: int = 5,
    max_pop_af: float = 1e-05,
    max_brca1_func_assay_score: float = -1.328,
    batch_size: int = 1000000,
):
    """
    Create all necessary views and tables in the DuckDB database.

    Sets up the database structure by creating views for filtering variants based on
    quality criteria, extracting information from the vals_info table, and preparing
    data for export.

    :param sample_id: the sample ID
    :param db: DuckDB database connection
    :param min_af: Minimum allele frequency threshold for variant filtering
    :param min_depth: Minimum read depth threshold for variant filtering
    :param max_pop_af: Maximum population allele frequency for variant filtering
    :param max_brca1_func_assay_score: Maximum BRCA1 functional assay score for variants
    :param batch_size: number of variants to process at a time
    """

    # delete lowest quality variants
    delete_low_quality_variants(db, min_af, min_depth)

    # create empty tables for batch processing
    make_supporting_tables(db)

    # make views needed to construct variants_enriched
    make_batch_views(db, limit=0, offset=0)
    make_vals_wide(db)
    make_info_wide(db)
    make_transcript_likely_lof(db)
    make_hgnc(db)
    make_oncogene_tsg(db)
    make_rescues(db, max_brca1_func_assay_score)

    # split the variants into evenly sized batches
    n_variants = db.table("variants").shape[0]
    n_batches = 1 + n_variants // batch_size
    chosen_batch_size = ceil(n_variants / n_batches)

    for i in range(n_batches):
        logging.info(f"Loading batch {i + 1} of {n_batches}")
        offset = i * chosen_batch_size

        # recreate views for this batch's variants, vals, and info
        make_batch_views(db, limit=chosen_batch_size, offset=offset)

        # populate fresh vep table, which is a table for performance reasons
        populate_vep(db)

        # append batch to variants_enriched
        make_variants_enriched(db, sample_id, max_pop_af)

    # populate mut_sig and somatic_variants from entire variants_enriched table
    make_mut_sig(db)
    make_somatic_variants(db)


def delete_low_quality_variants(
    db: duckdb.DuckDBPyConnection, min_af: float, min_depth: int
) -> None:
    """
    Delete variants based on allele frequency, read depth, and collected filters.

    :param db: DuckDB database connection
    :param min_af: Minimum allele frequency threshold
    :param min_depth: Minimum read depth threshold
    """

    logging.info("Deleting low quality variants")

    low_quality_vids = db.sql(f"""
        SELECT DISTINCT
            vid
        FROM
            vals_info
        WHERE (
            kind = 'val'
            AND k = 'af'
            AND v_float < {min_af}
        )
        OR (
            kind = 'val'
            AND k = 'dp'
            AND v_integer < {min_depth}
        )
        OR vid IN (
            SELECT
                vid
            FROM (
                SELECT
                    vid,
                    unnest(
                        str_split_regex(lower(v_varchar), '\s*\|\s*')
                    ) AS filter_val
                FROM
                    vals_info
                WHERE
                    kind = 'info'
                    AND k = 'as_filter_status'
                UNION
                SELECT
                    vid,
                    lower(unnest(filters)) AS filter_val
                FROM
                    variants
            )
            WHERE filter_val IN (
                'multiallelic',
                'map_qual',
                'slippage',
                'strand_bias',
                'weak_evidence', 
                'base_qual'
            )
        )
    """).fetchdf()

    db.register("low_quality_vids", low_quality_vids)

    db.sql("""
        DELETE FROM vals_info WHERE vid in (SELECT vid FROM low_quality_vids);
        DELETE FROM variants WHERE vid in (SELECT vid FROM low_quality_vids);
    """)

    db.unregister("low_quality_vids")


def make_batch_views(
    db: duckdb.DuckDBPyConnection, limit: int = 0, offset: int = 0
) -> None:
    """
    Create views that separate the vals_info table into 'vals' and 'info' views.

    :param db: DuckDB database connection
    :param limit: limit for select statement
    :param offset: offset for select statement
    """

    logging.info("Making vals and info views")

    db.sql(f"""
        CREATE OR REPLACE VIEW variants_batch AS (
            SELECT
                *
            FROM
                variants
            WHERE vid IN (
                SELECT vid
                FROM variants
                ORDER BY vid
                LIMIT {limit} OFFSET {offset}
            )
        )
    """)

    db.sql(f"""
        CREATE OR REPLACE VIEW vals AS (
            SELECT
                *
            FROM
                vals_info
            WHERE
                kind = 'val'
                AND vid IN (
                    SELECT vid
                    FROM variants
                    ORDER BY vid
                    LIMIT {limit} OFFSET {offset}
                )
        )
    """)

    db.sql(f"""
        CREATE OR REPLACE VIEW info AS (
            SELECT
                *
            FROM
                vals_info
            WHERE
                kind = 'info'
                AND vid IN (
                    SELECT vid
                    FROM variants
                    ORDER BY vid
                    LIMIT {limit} OFFSET {offset}
                )
        )
    """)


def make_supporting_tables(db: duckdb.DuckDBPyConnection) -> None:
    """
    Create empty physical tables to be populated later.

    :param db: DuckDB database connection
    """

    db.sql("""
        DROP TABLE IF EXISTS vep;
        
        CREATE TABLE vep (
            vid VARCHAR PRIMARY KEY,
            am_class VARCHAR,
            am_pathogenicity FLOAT,
            biotype VARCHAR,
            clin_sig VARCHAR,
            consequence VARCHAR,
            ensp VARCHAR,
            existing_variation VARCHAR,
            exon VARCHAR,
            feature VARCHAR,
            gene VARCHAR,
            gnom_ade_af FLOAT,
            gnom_adg_af FLOAT,
            hgnc_id VARCHAR,
            hgvsc VARCHAR,
            hgvsp VARCHAR,
            impact VARCHAR,
            intron VARCHAR,
            loftool FLOAT,
            mane_select VARCHAR,
            pli_gene_value FLOAT,
            poly_phen VARCHAR,
            sift VARCHAR,
            somatic VARCHAR,
            swissprot VARCHAR,
            symbol VARCHAR,
            uniprot_isoform VARCHAR,
            variant_class VARCHAR
        );
        
        DROP TABLE IF EXISTS variants_enriched;
        
        CREATE TABLE variants_enriched (
            sample_id VARCHAR,
            vid VARCHAR PRIMARY KEY,
            chrom VARCHAR,
            pos UINTEGER,
            id VARCHAR,
            ref VARCHAR,
            alt VARCHAR,
            qual VARCHAR,
            ref_count INTEGER,
            alt_count INTEGER,
            af FLOAT,
            dp INTEGER,
            gt VARCHAR,
            ps INTEGER,
            somatic BOOLEAN,
            impactful_splice_event BOOLEAN,
            in_clustered_event BOOLEAN,
            segdup BOOLEAN,
            repeat_masker BOOLEAN,
            low_pop_prevalence BOOLEAN,
            pon BOOLEAN,
            rescued BOOLEAN,
            rescued_cmc_tier BOOLEAN,
            rescued_hess BOOLEAN,
            rescued_met BOOLEAN,
            rescued_oc_brca1_func_assay_score BOOLEAN,
            rescued_oncogene_high_impact BOOLEAN,
            rescued_oncokb_hotspot BOOLEAN,
            rescued_oncokb_muteff BOOLEAN,
            rescued_oncokb_oncogenic BOOLEAN,
            rescued_tert BOOLEAN,
            rescued_tumor_suppressor_high_impact BOOLEAN,
            civic_desc VARCHAR,
            civic_id INTEGER,
            civic_score FLOAT,
            cosmic_tier INTEGER,
            gc_prop FLOAT,
            hess_driver BOOLEAN,
            hess_signature VARCHAR,
            lof JSON[],
            mc VARCHAR,
            nmd JSON[],
            oc_brca1_func_assay_score FLOAT,
            oc_gtex_gtex_gene VARCHAR,
            oc_gwas_catalog_disease VARCHAR,
            oc_gwas_catalog_pmid VARCHAR,
            oc_pharmgkb_id VARCHAR,
            oc_provean_prediction VARCHAR,
            oc_revel_score FLOAT,
            oncogene_high_impact BOOLEAN,
            oncokb_hotspot BOOLEAN,
            oncokb_muteff VARCHAR,
            oncokb_oncogenic VARCHAR,
            protein_changed BOOLEAN,
            rs VARCHAR,
            tumor_suppressor_high_impact BOOLEAN,
            vep_am_class VARCHAR,
            vep_am_pathogenicity FLOAT,
            vep_biotype VARCHAR,
            vep_clin_sig VARCHAR,
            vep_consequence VARCHAR,
            vep_ensp VARCHAR,
            vep_existing_variation VARCHAR,
            vep_exon VARCHAR,
            vep_feature VARCHAR,
            vep_gene VARCHAR,
            vep_gnom_ade_af FLOAT,
            vep_gnom_adg_af FLOAT,
            vep_hgnc_id VARCHAR,
            vep_hgvsc VARCHAR,
            vep_hgvsp VARCHAR,
            vep_impact VARCHAR,
            vep_intron VARCHAR,
            vep_loftool FLOAT,
            vep_mane_select VARCHAR,
            vep_pli_gene_value FLOAT,
            vep_poly_phen VARCHAR,
            vep_sift VARCHAR,
            vep_somatic VARCHAR,
            vep_swissprot VARCHAR,
            vep_symbol VARCHAR,
            vep_uniprot_isoform VARCHAR,
            vep_variant_class VARCHAR
        );
        
        DROP TABLE IF EXISTS somatic_variants;
        
        CREATE TABLE somatic_variants (
            sample_id VARCHAR,
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
            lof VARCHAR,
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


def make_vals_wide(db: duckdb.DuckDBPyConnection) -> None:
    """
    Create a wide view of data from the vals view.

    Transforms the vals view from a long format to a wide format, with columns for
    ref/alt counts, allele frequency, read depth, genotype, and phase set.

    :param db: DuckDB database connection
    """

    logging.info("Making vals_wide")

    db.sql("""
        CREATE OR REPLACE VIEW vals_wide AS (
            SELECT
                DISTINCT vals.vid,
                ref_count: t_ad.ad[1],
                alt_count: t_ad.ad[2],
                t_af.af,
                t_dp.dp,
                t_gt.gt,
                t_ps.ps
            FROM
                vals

            LEFT JOIN (
                SELECT
                    vid,
                    ad: v_integer_arr
                FROM
                    vals
                WHERE
                    k = 'ad'
            ) t_ad ON vals.vid = t_ad.vid

            LEFT JOIN (
                SELECT
                    vid,
                    af: v_float
                FROM
                    vals
                WHERE
                    k = 'af'
            ) t_af ON vals.vid = t_af.vid

            LEFT JOIN (
                SELECT
                    vid,
                    dp: v_integer
                FROM
                    vals
                WHERE
                    k = 'dp'
            ) t_dp ON vals.vid = t_dp.vid

            LEFT JOIN (
                SELECT
                    vid,
                    gt: v_varchar
                FROM
                    vals
                WHERE
                    k = 'gt'
            ) t_gt ON vals.vid = t_gt.vid

            LEFT JOIN (
                SELECT
                    vid,
                    ps: v_integer
                FROM
                    vals
                WHERE
                    k = 'ps'
            ) t_ps ON vals.vid = t_ps.vid
        )
    """)


def make_info_wide(db: duckdb.DuckDBPyConnection) -> None:
    """
    Create a wide view of variant information from the info table.

    Transforms the info table from a long format to a wide format, with columns for
    various annotations like BRCA1 functional assay score, CIViC information,
    COSMIC tier, dbSNP ID, and other annotations.

    :param db: DuckDB database connection
    """

    logging.info("Making info_wide")

    db.sql("""
        CREATE OR REPLACE VIEW info_wide AS (
            SELECT
                DISTINCT info.vid,
                t_oc_brca1_func_assay_score.oc_brca1_func_assay_score,
                t_civic_desc.civic_desc,
                t_civic_id.civic_id,
                t_civic_score.civic_score,
                t_cosmic_tier.cosmic_tier,
                rs: t_rs.rs[1],
                t_gc_prop.gc_prop,
                mc: list_aggregate(t_mc.mc, 'string_agg', ','),
                t_lof.lof,
                t_nmd.nmd,
                t_oc_gtex_gtex_gene.oc_gtex_gtex_gene,
                t_oc_gwas_catalog_disease.oc_gwas_catalog_disease,
                t_oc_gwas_catalog_pmid.oc_gwas_catalog_pmid,
                hess_driver: coalesce(t_hess.hess_driver, FALSE),
                t_hess.hess_signature,
                segdup: coalesce(t_segdup.segdup, FALSE),
                repeat_masker: coalesce(t_repeat_masker.repeat_masker, FALSE),
                pon: coalesce(t_pon.pon, FALSE),
                t_oc_pharmgkb_id.oc_pharmgkb_id,
                t_oc_provean_prediction.oc_provean_prediction,
                t_oc_revel_score.oc_revel_score,
                oncokb_hotspot: coalesce(t_oncokb_hotspot.oncokb_hotspot, FALSE),
                t_oncokb_muteff.oncokb_muteff,
                t_oncokb_oncogenic.oncokb_oncogenic
            FROM
                info

            LEFT JOIN (
                SELECT
                    vid,
                    oc_brca1_func_assay_score: v_float
                FROM
                    info
                WHERE
                    k = 'oc_brca1_func_assay_score'
            ) t_oc_brca1_func_assay_score ON
                info.vid = t_oc_brca1_func_assay_score.vid

            LEFT JOIN (
                SELECT
                    vid,
                    civic_desc: v_varchar
                FROM
                    info
                WHERE
                    k = 'civic_desc'
            ) t_civic_desc ON info.vid = t_civic_desc.vid

            LEFT JOIN (
                SELECT
                    vid,
                    civic_id: v_integer
                FROM
                    info
                WHERE
                    k = 'civic_id'
            ) t_civic_id ON info.vid = t_civic_id.vid

            LEFT JOIN (
                SELECT
                    vid,
                    civic_score: v_float
                FROM
                    info
                WHERE
                    k = 'civic_score'
            ) t_civic_score ON info.vid = t_civic_score.vid

            LEFT JOIN (
                SELECT
                    vid,
                    cosmic_tier: v_integer
                FROM
                    info
                WHERE
                    k = 'cmc_tier'
            ) t_cosmic_tier ON info.vid = t_cosmic_tier.vid

            LEFT JOIN (
                SELECT
                    vid,
                    rs: v_varchar_arr
                FROM
                    info
                WHERE
                    k = 'rs'
            ) t_rs ON info.vid = t_rs.vid

            LEFT JOIN (
                SELECT
                    vid,
                    gc_prop: v_float
                FROM
                    info
                WHERE
                    k = 'gc_prop'
            ) t_gc_prop ON info.vid = t_gc_prop.vid

            LEFT JOIN (
                SELECT
                    vid,
                    mc: v_varchar_arr
                FROM
                    info
                WHERE
                    k = 'mc'
            ) t_mc ON info.vid = t_mc.vid

            LEFT JOIN (
                SELECT
                    vid,
                    lof: v_json_arr
                FROM
                    info
                WHERE
                    k = 'lof'
            ) t_lof ON info.vid = t_lof.vid

            LEFT JOIN (
                SELECT
                    vid,
                    nmd: v_json_arr
                FROM
                    info
                WHERE
                    k = 'nmd'
            ) t_nmd ON info.vid = t_nmd.vid

            LEFT JOIN (
                SELECT
                    vid,
                    oc_gtex_gtex_gene: v_varchar
                FROM
                    info
                WHERE
                    k = 'oc_gtex_gtex_gene'
            ) t_oc_gtex_gtex_gene ON info.vid = t_oc_gtex_gtex_gene.vid

            LEFT JOIN (
                SELECT
                    vid,
                    oc_gwas_catalog_disease: v_varchar
                FROM
                    info
                WHERE
                    k = 'oc_gwas_catalog_disease'
            ) t_oc_gwas_catalog_disease ON info.vid = t_oc_gwas_catalog_disease.vid

            LEFT JOIN (
                SELECT
                    vid,
                    oc_gwas_catalog_pmid: v_varchar
                FROM
                    info
                WHERE
                    k = 'oc_gwas_catalog_pmid'
            ) t_oc_gwas_catalog_pmid ON info.vid = t_oc_gwas_catalog_pmid.vid

            LEFT JOIN (
                SELECT
                    vid,
                    hess_signature: v_varchar,
                    hess_driver: TRUE
                FROM
                    info
                WHERE
                    k = 'hess'
            ) t_hess ON info.vid = t_hess.vid

            LEFT JOIN (
                SELECT
                    vid,
                    segdup: TRUE
                FROM
                    info
                WHERE
                    k = 'segdup'
            ) t_segdup ON info.vid = t_segdup.vid

            LEFT JOIN (
                SELECT
                    vid,
                    repeat_masker: TRUE
                FROM
                    info
                WHERE
                    k = 'rm'
            ) t_repeat_masker ON info.vid = t_repeat_masker.vid

            LEFT JOIN (
                SELECT
                    vid,
                    pon: TRUE
                FROM
                    info
                WHERE
                    k = 'pon'
            ) t_pon ON info.vid = t_pon.vid

            LEFT JOIN (
                SELECT
                    vid,
                    oc_pharmgkb_id: v_varchar
                FROM
                    info
                WHERE
                    k = 'oc_pharmgkb_id'
            ) t_oc_pharmgkb_id ON info.vid = t_oc_pharmgkb_id.vid

            LEFT JOIN (
                SELECT
                    vid,
                    oc_provean_prediction: v_varchar
                FROM
                    info
                WHERE
                    k = 'oc_provean_prediction'
            ) t_oc_provean_prediction ON info.vid = t_oc_provean_prediction.vid

            LEFT JOIN (
                SELECT
                    vid,
                    oc_revel_score: v_float
                FROM
                    info
                WHERE
                    k = 'oc_revel_score'
            ) t_oc_revel_score ON info.vid = t_oc_revel_score.vid

            LEFT JOIN (
                SELECT
                    vid,
                    oncokb_muteff: v_varchar
                FROM
                    info
                WHERE
                    k = 'oncokb_muteff'
            ) t_oncokb_muteff ON info.vid = t_oncokb_muteff.vid

            LEFT JOIN (
                SELECT
                    vid,
                    oncokb_hotspot: v_boolean
                FROM
                    info
                WHERE
                    k = 'oncokb_hotspot'
            ) t_oncokb_hotspot ON info.vid = t_oncokb_hotspot.vid

            LEFT JOIN (
                SELECT
                    vid,
                    oncokb_oncogenic: v_varchar
                FROM
                    info
                WHERE
                    k = 'oncokb_oncogenic'
            ) t_oncokb_oncogenic ON info.vid = t_oncokb_oncogenic.vid
        )
    """)


def make_transcript_likely_lof(db: duckdb.DuckDBPyConnection) -> None:
    """
    Create a view of transcripts likely to have loss-of-function (LoF) variants.

    Processes the oc_revel_all field from the info table to identify transcripts with a
    REVEL score >= 0.7, which indicates likely loss of function.

    :param db: DuckDB database connection
    """

    logging.info("Making transcript_likely_lof_v")

    db.sql("""
        CREATE OR REPLACE VIEW transcript_likely_lof_v AS (
            WITH exploded AS (
                SELECT
                    vid,
                    revel_cols: unnest(
                        json_transform_strict(v_varchar, '["json"]')
                    )
                FROM
                    info
                WHERE
                    k = 'oc_revel_all'
            ),
            split_to_cols AS (
                SELECT
                    vid,
                    transcript_id: json_extract_string(revel_cols, '$[0]'),
                    score: json_extract(revel_cols, '$[1]')::DOUBLE
                FROM exploded
            )
            SELECT
                vid,
                transcript_likely_lof: string_agg(transcript_id, ';')
            FROM
                split_to_cols
            WHERE
                score >= 0.7
            GROUP BY
                vid
        )
    """)


def make_hgnc(db: duckdb.DuckDBPyConnection) -> None:
    """
    Create a view of HGNC identifiers.

    :param db: DuckDB database connection
    """

    logging.info("Making hgnc_v")

    db.sql("""
        CREATE OR REPLACE VIEW hgnc_v AS (
            WITH hgnc_name_exploded AS (
                SELECT
                    vid,
                    v_varchar: unnest(v_varchar_arr)
                FROM
                    info
                WHERE
                    k = 'hgnc_name'
            ),
            hgnc_name_concat AS (
                SELECT
                    vid,
                    hgnc_name: string_agg(v_varchar, ';')
                FROM
                    hgnc_name_exploded
                GROUP BY
                    vid
            ),
            hgnc_group_exploded AS (
                SELECT
                    vid,
                    v_varchar: unnest(v_varchar_arr)
                FROM
                    info
                WHERE
                    k = 'hgnc_group'
            ),
            hgnc_group_concat AS (
                SELECT
                    vid,
                    hgnc_group: string_agg(v_varchar, ';')
                FROM
                    hgnc_group_exploded
                GROUP BY
                    vid
            )
            SELECT
                vid: coalesce(hgnc_name_concat.vid, hgnc_group_concat.vid),
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


def populate_vep(db: duckdb.DuckDBPyConnection) -> None:
    """
    Create a table of Variant Effect Predictor (VEP) annotations.

    Processes the csq field from the info table to extract VEP annotations for each
    variant. This is a physical table due to the relative slowness of the
    `vep_exploded` CTE, the fact it gets queried several times later, and that views in
    DuckDB aren't materialized.

    :param db: DuckDB database connection
    """

    logging.info("Populating vep")

    db.sql("TRUNCATE vep;")

    db.sql("""
        INSERT INTO vep BY NAME (
            WITH vep_exploded AS (
                SELECT
                    vid,
                    csq: json_transform(
                        unnest(v_json_arr),
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
                    )
                FROM
                    info
                WHERE
                    k = 'csq'
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
                hgvsp: url_decode(csq.hgvsp),
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


def make_oncogene_tsg(db: duckdb.DuckDBPyConnection) -> None:
    """
    Create a view of oncogene and tumor suppressor gene (TSG) information.

    Identifies variants in oncogenes and tumor suppressor genes with high-impact
    effects, based on the oncogene and tsg fields from the info table and the impact
    field from the vep table.

    :param db: DuckDB database connection
    """

    logging.info("Making oncogene_tsg")

    db.sql("""
        CREATE OR REPLACE VIEW oncogene_tsg AS (
            WITH oncogene_v AS (
                SELECT
                    vid,
                    oncogene_high_impact: TRUE
                FROM
                    info
                WHERE
                    k = 'oncogene'
                    AND vid IN (SELECT vid FROM vep WHERE impact = 'HIGH')
            ),
            tsg_v AS (
                SELECT
                    vid,
                    tumor_suppressor_high_impact: TRUE
                FROM
                    info
                WHERE
                    k = 'tsg'
                    AND vid IN (SELECT vid FROM vep WHERE impact = 'HIGH')
            )
            SELECT
                vid: coalesce(oncogene_v.vid, tsg_v.vid),
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


def make_rescues(
    db: duckdb.DuckDBPyConnection, max_brca1_func_assay_score: float
) -> None:
    """
    Create a view of rescued variants based on various criteria.

    Identifies variants that should be rescued despite not meeting standard quality
    criteria, based on factors like oncogenic effect, hotspot status, COSMIC tier,
    BRCA1 functional assay score, oncogene/TSG status, and specific gene-position
    combinations (e.g., TERT promoter, MET exon 14 skipping).

    :param db: DuckDB database connection
    :param max_brca1_func_assay_score: Maximum BRCA1 functional assay score for variants
    """

    logging.info("Making rescues")

    db.sql(f"""
        CREATE OR REPLACE VIEW rescues AS (
            WITH rescue_reasons AS (
                SELECT
                    vid,
                    rescued_oncokb_muteff: (
                        oncokb_muteff IN ('Loss-of-function', 'Gain-of-function')
                    ),
                    rescued_oncokb_oncogenic: oncokb_oncogenic = 'Oncogenic',
                    rescued_oncokb_hotspot: oncokb_hotspot,
                    rescued_cmc_tier: cosmic_tier = 1,
                    rescued_oc_brca1_func_assay_score: (
                        oc_brca1_func_assay_score <= {max_brca1_func_assay_score}
                    ),
                    rescued_oncogene_high_impact: vid IN (
                        SELECT
                            vid
                        FROM
                            oncogene_tsg
                        WHERE
                            oncogene_high_impact
                    ),
                    rescued_tumor_suppressor_high_impact: vid IN (
                        SELECT
                            vid
                        FROM
                            oncogene_tsg
                        WHERE
                            tumor_suppressor_high_impact
                    ),
                    rescued_hess: hess_driver,
                    rescued_tert: vid IN (
                        SELECT
                            variants_batch.vid
                        FROM
                            variants_batch
                        INNER JOIN
                            vep
                        ON
                            variants_batch.vid = vep.vid
                        WHERE
                            variants_batch.chrom = 'chr5'
                            AND
                            variants_batch.pos BETWEEN 1295054 AND 1295365
                            AND
                            vep.symbol = 'TERT'
                    ),
                    rescued_met: vid IN (
                        SELECT
                            variants_batch.vid
                        FROM
                            variants_batch
                        INNER JOIN
                            vep
                        ON
                            variants_batch.vid = vep.vid
                        WHERE
                            variants_batch.chrom = 'chr7'
                            AND
                            variants_batch.pos BETWEEN 116771825 AND 116771840
                            AND
                            vep.symbol = 'MET'
                    )
                FROM info_wide
            )
            SELECT
                *,
                rescued: TRUE
            FROM
                rescue_reasons
            WHERE
                rescued_oncokb_muteff
                OR rescued_oncokb_oncogenic
                OR rescued_oncokb_hotspot
                OR rescued_cmc_tier
                OR rescued_oc_brca1_func_assay_score
                OR rescued_oncogene_high_impact
                OR rescued_tumor_suppressor_high_impact
                OR rescued_hess
                OR rescued_tert
                OR rescued_met
        )
    """)


def make_variants_enriched(
    db: duckdb.DuckDBPyConnection, sample_id: str, max_pop_af: float
) -> None:
    """
    Create an enriched, wide view of variants, assigning a `somatic` column if the
    variant is identified as a high-quality somatic variant.

    Identifies variants that pass quality filters and either are rescued or are valid
    somatic alterations (e.g. splice events, protein changes) while not being in
    clustered events, segmental duplications, or repeat regions, and having low
    population allele frequency.

    :param sample_id: the sample ID
    :param db: DuckDB database connection
    :param max_pop_af: Maximum population allele frequency
    """

    logging.info("Populating variants_enriched")

    db.sql(f"""
        INSERT INTO variants_enriched BY NAME (
            WITH variants_wide AS (
                SELECT
                    variants_batch.* EXCLUDE filters,
                    vals_wide.* EXCLUDE vid,
                    info_wide.* EXCLUDE vid,
                    COLUMNS(vep.*) AS "vep_\\0",
                    oncogene_high_impact: coalesce(
                        oncogene_tsg.oncogene_high_impact, FALSE
                    ),
                    tumor_suppressor_high_impact: coalesce(
                        oncogene_tsg.tumor_suppressor_high_impact, FALSE
                    ),
                    rescued_oncokb_muteff: coalesce(
                        rescues.rescued_oncokb_muteff, FALSE
                    ),
                    rescued_oncokb_oncogenic: coalesce(
                        rescues.rescued_oncokb_oncogenic, FALSE
                    ),
                    rescued_oncokb_hotspot: coalesce(
                        rescues.rescued_oncokb_hotspot, FALSE
                    ),
                    rescued_cmc_tier: coalesce(rescues.rescued_cmc_tier, FALSE),
                    rescued_oc_brca1_func_assay_score: coalesce(
                        rescues.rescued_oc_brca1_func_assay_score, FALSE
                    ),
                    rescued_oncogene_high_impact: coalesce(
                        rescues.rescued_oncogene_high_impact, FALSE
                    ),
                    rescued_tumor_suppressor_high_impact: coalesce(
                        rescues.rescued_tumor_suppressor_high_impact, FALSE
                    ),
                    rescued_hess: coalesce(rescues.rescued_hess, FALSE),
                    rescued_tert: coalesce(rescues.rescued_tert, FALSE),
                    rescued_met: coalesce(rescues.rescued_met, FALSE),
                    rescued: coalesce(rescues.rescued, FALSE)
                FROM
                    variants_batch
                    
                LEFT JOIN vals_wide ON variants_batch.vid = vals_wide.vid
                LEFT JOIN info_wide ON variants_batch.vid = info_wide.vid
                LEFT JOIN hgnc_v ON variants_batch.vid = hgnc_v.vid
                LEFT JOIN vep ON variants_batch.vid = vep.vid
                LEFT JOIN oncogene_tsg ON variants_batch.vid = oncogene_tsg.vid
                LEFT JOIN rescues ON variants_batch.vid = rescues.vid
            ),
            summarized AS (
                SELECT
                    * EXCLUDE vep_vid,
                    impactful_splice_event: (
                        contains(vep_consequence, 'splice') AND
                        vep_impact IN ('HIGH', 'MODERATE')
                    ),
                    protein_changed: (
                        vep_hgvsp IS NOT NULL
                        AND
                        NOT ends_with(vep_hgvsp, '=')
                    ),
                    in_clustered_event: vid IN (
                        SELECT
                            vid
                        FROM (
                            SELECT
                                vid,
                                lower(unnest(filters)) AS filter_val
                            FROM
                            variants_batch
                        )
                        WHERE
                            filter_val = 'clustered_events'
                    ),
                    low_pop_prevalence: (
                        coalesce(vep_gnom_ade_af, 0) <= {max_pop_af}
                        AND
                        coalesce(vep_gnom_adg_af, 0) <= {max_pop_af}
                        AND
                        NOT pon
                    )
                FROM variants_wide
            )
            SELECT
                sample_id: '{sample_id}',
                *,
                somatic: (
                    rescued
                    OR (
                        (impactful_splice_event OR protein_changed)
                        AND NOT (in_clustered_event OR segdup OR repeat_masker)
                        AND low_pop_prevalence
                    )
                )
            FROM
                summarized
        )
    """)


def make_mut_sig(db: duckdb.DuckDBPyConnection) -> None:
    """
    Create a minimal view of enriched_variants for use in mutational signature analysis.

    :param db: DuckDB database connection
    """

    logging.info("Creating mut_sig")

    db.sql(f"""
        CREATE OR REPLACE VIEW mut_sig AS (
            SELECT
                sample_id,
                chrom,
                pos,
                id,
                ref,
                alt,
                ref_count,
                alt_count,
                af,
                dp,
                rs,
                segdup: coalesce(segdup, FALSE),
                repeat_masker: coalesce(repeat_masker, FALSE),
                pon: coalesce(pon, FALSE),
                in_clustered_event: coalesce(in_clustered_event, FALSE),
                cosmic_tier,
                oncokb_oncogenic,
                protein_changed,
                vep_biotype,
                vep_consequence,
                vep_gnom_ade_af,
                vep_gnom_adg_af,
                vep_variant_class
            FROM
                variants_enriched
            ORDER BY
                replace(replace(chrom[4:], 'X', '23'), 'Y', '24')::INTEGER,
                pos,
                alt
        )
    """)


def make_somatic_variants(db: duckdb.DuckDBPyConnection) -> None:
    """
    Create a table of somatic variants with all annotations.

    Combines information from all previously created views and tables to create a
    comprehensive table of somatic variants with all annotations.

    :param db: DuckDB database connection
    """

    logging.info("Populating somatic_variants")

    db.sql(f"""
        INSERT INTO
            somatic_variants
        BY NAME (
            SELECT
                sample_id,
                chrom,
                pos,
                ref,
                alt,
                ref_count,
                alt_count,
                af,
                dp,
                gt,
                ps,
                brca1_func_score: oc_brca1_func_assay_score,
                civic_description: civic_desc,
                civic_id,
                civic_score,
                cosmic_tier,
                dbsnp_rs_id: rs,
                gc_content: gc_prop,
                lof,
                molecular_consequence: mc,
                nmd,
                gtex_gene: oc_gtex_gtex_gene,
                gwas_disease: oc_gwas_catalog_disease,
                gwas_pmid: oc_gwas_catalog_pmid,
                pharmgkb_id: oc_pharmgkb_id,
                provean_prediction: oc_provean_prediction,
                revel_score: oc_revel_score,
                oncokb_effect: oncokb_muteff,
                oncokb_hotspot,
                oncokb_oncogenic,
                hess_driver,
                hess_signature,
                am_class: vep_am_class,
                am_pathogenicity: vep_am_pathogenicity,
                dna_change: vep_hgvsc,
                ensembl_feature_id: vep_feature,
                ensembl_gene_id: vep_gene,
                exon: vep_exon,
                gnomade_af: vep_gnom_ade_af,
                gnomadg_af: vep_gnom_adg_af,
                hugo_symbol: vep_symbol,
                intron: vep_intron,
                polyphen: vep_poly_phen,
                protein_change: vep_hgvsp,
                sift: vep_sift,
                uniprot_id: vep_uniprot_isoform,
                variant_info: vep_consequence,
                variant_type: vep_variant_class,
                vep_biotype,
                vep_clin_sig,
                vep_ensp,
                vep_existing_variation,
                vep_hgnc_id,
                vep_impact,
                vep_loftool,
                vep_mane_select,
                vep_pli_gene_value,
                vep_somatic,
                vep_swissprot,
                transcript_likely_lof_v.transcript_likely_lof,
                hgnc_v.hgnc_name,
                hgnc_family: hgnc_v.hgnc_group,
                oncogene_high_impact,
                tumor_suppressor_high_impact,
                rescue: rescued
            FROM
                variants_enriched
            
            LEFT JOIN
                transcript_likely_lof_v
            ON
                variants_enriched.vid = transcript_likely_lof_v.vid
            
            LEFT JOIN
                hgnc_v
            ON
                variants_enriched.vid = hgnc_v.vid
            
            WHERE somatic
            
            ORDER BY
                replace(replace(chrom[4:], 'X', '23'), 'Y', '24')::INTEGER,
                pos,
                alt
        )
    """)


def get_somatic_variants_as_df(db: duckdb.DuckDBPyConnection) -> pd.DataFrame:
    """
    Retrieve somatic variants from the database as a pandas DataFrame.

    :param db: DuckDB database connection
    :returns: A pandas DataFrame containing somatic variants with all annotations
    """

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

    return (
        db.table("somatic_variants")
        .df()
        .astype(
            {
                "sample_id": "string",
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
                "lof": "string",
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
