from tempfile import NamedTemporaryFile, TemporaryDirectory
from typing import Generator

import duckdb
import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

from annotate_mutations_postprocess.maf import get_somatic_variants_as_df, make_views

maf_dtypes = {
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


@pytest.fixture(scope="class")
def db_setup() -> Generator:
    with TemporaryDirectory() as tmpdir:
        db_path = NamedTemporaryFile(dir=tmpdir, suffix="duckdb").name

        with duckdb.connect(db_path) as db:
            db.sql("""
                CREATE TABLE variants(
                    vid UINTEGER PRIMARY KEY,
                    chrom VARCHAR NOT NULL,
                    pos UINTEGER NOT NULL,
                    id VARCHAR,
                    "ref" VARCHAR,
                    alt VARCHAR,
                    qual VARCHAR,
                    filters VARCHAR[]
                );
                
                CREATE TABLE vals_info(
                    vid UINTEGER REFERENCES variants (vid),
                    kind VARCHAR,
                    k VARCHAR NOT NULL,
                    v_boolean BOOLEAN,
                    v_varchar VARCHAR,
                    v_integer INTEGER,
                    v_float FLOAT,
                    v_json JSON,
                    v_boolean_arr BOOLEAN[],
                    v_varchar_arr VARCHAR[],
                    v_integer_arr INTEGER[],
                    v_float_arr FLOAT[],
                    v_json_arr JSON[]
                );
            """)

            yield db


@pytest.fixture
def db(db_setup: duckdb.DuckDBPyConnection) -> Generator:
    db_setup.sql("""
        DROP TABLE IF EXISTS somatic_variants;
        DROP TABLE IF EXISTS vep;
        
        TRUNCATE TABLE vals_info;
        TRUNCATE TABLE variants;
    """)

    yield db_setup


@pytest.mark.usefixtures("db_setup")
class Test:
    def test_noop(self, db: duckdb.DuckDBPyConnection) -> None:
        make_views(db)

        observed = get_somatic_variants_as_df(db)
        expected = pd.DataFrame(
            [],
            columns=list(maf_dtypes.keys()),  # pyright: ignore
        ).astype(maf_dtypes)
        assert_frame_equal(observed, expected)


@pytest.mark.usefixtures("db_setup")
class TestFilter:
    def test_filter(self, db: duckdb.DuckDBPyConnection) -> None:
        db.sql("""
            INSERT INTO variants(
                vid, chrom, pos, ref, alt, filters
            ) VALUES (
                1, 'chr1', 100, 'G', 'C', ['germline']
            ), (
                2, 'chr2', 200, 'A', 'T', ['slippage', 'irrelevant']
            ), (
                3, 'chr3', 300, 'A', 'T', ['PASS']
            );
            
            INSERT INTO vals_info(
                vid, kind, k, v_integer, v_float
            ) VALUES (
                1, 'val', 'af', NULL, 0.3
            ), (
                1, 'val', 'dp', 20, NULL
            ), (
                2, 'val', 'af', NULL, 0.3
            ), (
                2, 'val', 'dp', 20, NULL
            ), (
                3, 'val', 'af', NULL, 0.3
            ), (
                3, 'val', 'dp', 20, NULL
            );
            
            --any rescue state will work for this test
            INSERT INTO vals_info(
                vid, kind, k, v_varchar
            ) VALUES (
                1, 'info', 'oncokb_muteff', 'Loss-of-function'
            ), (
                1, 'info', 'as_filter_status', 'SITE|SITE'
            ), (
                2, 'info', 'oncokb_muteff', 'Loss-of-function'
            ), (
                3, 'info', 'oncokb_muteff', 'Loss-of-function'
            ), (
                3,'info', 'as_filter_status', 'strand_bias|SITE'
            );
        """)

        make_views(db)

        observed = get_somatic_variants_as_df(db)

        expected = pd.DataFrame(
            [
                {
                    "chrom": "chr1",
                    "pos": 100,
                    "ref": "G",
                    "alt": "C",
                    "af": 0.3,
                    "alt_count": None,
                    "am_class": None,
                    "am_pathogenicity": None,
                    "brca1_func_score": None,
                    "civic_description": None,
                    "civic_id": None,
                    "civic_score": None,
                    "cosmic_tier": None,
                    "dbsnp_rs_id": None,
                    "dna_change": None,
                    "dp": 20,
                    "ensembl_feature_id": None,
                    "ensembl_gene_id": None,
                    "exon": None,
                    "gc_content": None,
                    "gnomade_af": None,
                    "gnomadg_af": None,
                    "gt": None,
                    "gtex_gene": None,
                    "gwas_disease": None,
                    "gwas_pmid": None,
                    "hess_driver": False,
                    "hess_signature": None,
                    "hgnc_family": None,
                    "hgnc_name": None,
                    "hugo_symbol": None,
                    "intron": None,
                    "lof": None,
                    "molecular_consequence": None,
                    "nmd": None,
                    "oncogene_high_impact": False,
                    "oncokb_effect": "Loss-of-function",
                    "oncokb_hotspot": False,
                    "oncokb_oncogenic": None,
                    "pharmgkb_id": None,
                    "polyphen": None,
                    "protein_change": None,
                    "provean_prediction": None,
                    "ps": None,
                    "ref_count": None,
                    "rescue": True,
                    "revel_score": None,
                    "sift": None,
                    "transcript_likely_lof": None,
                    "tumor_suppressor_high_impact": False,
                    "uniprot_id": None,
                    "variant_info": None,
                    "variant_type": None,
                    "vep_biotype": None,
                    "vep_clin_sig": None,
                    "vep_ensp": None,
                    "vep_existing_variation": None,
                    "vep_hgnc_id": None,
                    "vep_impact": None,
                    "vep_loftool": None,
                    "vep_mane_select": None,
                    "vep_pli_gene_value": None,
                    "vep_somatic": None,
                    "vep_swissprot": None,
                },
            ]
        ).astype(maf_dtypes)

        assert_frame_equal(observed, expected)


@pytest.mark.usefixtures("db_setup")
class TestRescue:
    def test_oncokb_muteff(self, db: duckdb.DuckDBPyConnection) -> None:
        db.sql("""
            INSERT INTO variants(
                vid, chrom, pos, ref, alt
            ) VALUES (
                1, 'chr1', 100, 'G', 'C'
            ), (
                2, 'chr2', 200, 'A', 'T'
            ), (
                3, 'chr3', 300, 'A', 'T'
            );
            
            INSERT INTO vals_info(
                vid, kind, k, v_integer, v_float
            ) VALUES (
                1, 'val', 'af', NULL, 0.3
            ), (
                1, 'val', 'dp', 20, NULL
            ), (
                2, 'val', 'af', NULL, 0.3
            ), (
                2, 'val', 'dp', 20, NULL
            ), (
                3, 'val', 'af', NULL, 0.3
            ), (
                3, 'val', 'dp', 20, NULL
            );
            
            INSERT INTO vals_info(
                vid, kind, k, v_varchar
            ) VALUES (
                1, 'info', 'oncokb_muteff', 'Loss-of-function'
            ), (
                2, 'info', 'oncokb_muteff', 'Gain-of-function'
            ), (
                3, 'info', 'oncokb_muteff', 'other'
            );
        """)

        make_views(db)

        observed = get_somatic_variants_as_df(db)

        expected = pd.DataFrame(
            [
                {
                    "chrom": "chr1",
                    "pos": 100,
                    "ref": "G",
                    "alt": "C",
                    "af": 0.3,
                    "alt_count": None,
                    "am_class": None,
                    "am_pathogenicity": None,
                    "brca1_func_score": None,
                    "civic_description": None,
                    "civic_id": None,
                    "civic_score": None,
                    "cosmic_tier": None,
                    "dbsnp_rs_id": None,
                    "dna_change": None,
                    "dp": 20,
                    "ensembl_feature_id": None,
                    "ensembl_gene_id": None,
                    "exon": None,
                    "gc_content": None,
                    "gnomade_af": None,
                    "gnomadg_af": None,
                    "gt": None,
                    "gtex_gene": None,
                    "gwas_disease": None,
                    "gwas_pmid": None,
                    "hess_driver": False,
                    "hess_signature": None,
                    "hgnc_family": None,
                    "hgnc_name": None,
                    "hugo_symbol": None,
                    "intron": None,
                    "lof": None,
                    "molecular_consequence": None,
                    "nmd": None,
                    "oncogene_high_impact": False,
                    "oncokb_effect": "Loss-of-function",
                    "oncokb_hotspot": False,
                    "oncokb_oncogenic": None,
                    "pharmgkb_id": None,
                    "polyphen": None,
                    "protein_change": None,
                    "provean_prediction": None,
                    "ps": None,
                    "ref_count": None,
                    "rescue": True,
                    "revel_score": None,
                    "sift": None,
                    "transcript_likely_lof": None,
                    "tumor_suppressor_high_impact": False,
                    "uniprot_id": None,
                    "variant_info": None,
                    "variant_type": None,
                    "vep_biotype": None,
                    "vep_clin_sig": None,
                    "vep_ensp": None,
                    "vep_existing_variation": None,
                    "vep_hgnc_id": None,
                    "vep_impact": None,
                    "vep_loftool": None,
                    "vep_mane_select": None,
                    "vep_pli_gene_value": None,
                    "vep_somatic": None,
                    "vep_swissprot": None,
                },
                {
                    "chrom": "chr2",
                    "pos": 200,
                    "ref": "A",
                    "alt": "T",
                    "af": 0.3,
                    "alt_count": None,
                    "am_class": None,
                    "am_pathogenicity": None,
                    "brca1_func_score": None,
                    "civic_description": None,
                    "civic_id": None,
                    "civic_score": None,
                    "cosmic_tier": None,
                    "dbsnp_rs_id": None,
                    "dna_change": None,
                    "dp": 20,
                    "ensembl_feature_id": None,
                    "ensembl_gene_id": None,
                    "exon": None,
                    "gc_content": None,
                    "gnomade_af": None,
                    "gnomadg_af": None,
                    "gt": None,
                    "gtex_gene": None,
                    "gwas_disease": None,
                    "gwas_pmid": None,
                    "hess_driver": False,
                    "hess_signature": None,
                    "hgnc_family": None,
                    "hgnc_name": None,
                    "hugo_symbol": None,
                    "intron": None,
                    "lof": None,
                    "molecular_consequence": None,
                    "nmd": None,
                    "oncogene_high_impact": False,
                    "oncokb_effect": "Gain-of-function",
                    "oncokb_hotspot": False,
                    "oncokb_oncogenic": None,
                    "pharmgkb_id": None,
                    "polyphen": None,
                    "protein_change": None,
                    "provean_prediction": None,
                    "ps": None,
                    "ref_count": None,
                    "rescue": True,
                    "revel_score": None,
                    "sift": None,
                    "transcript_likely_lof": None,
                    "tumor_suppressor_high_impact": False,
                    "uniprot_id": None,
                    "variant_info": None,
                    "variant_type": None,
                    "vep_biotype": None,
                    "vep_clin_sig": None,
                    "vep_ensp": None,
                    "vep_existing_variation": None,
                    "vep_hgnc_id": None,
                    "vep_impact": None,
                    "vep_loftool": None,
                    "vep_mane_select": None,
                    "vep_pli_gene_value": None,
                    "vep_somatic": None,
                    "vep_swissprot": None,
                },
            ]
        ).astype(maf_dtypes)

        assert_frame_equal(observed, expected)

    def test_oncokb_oncogenic(self, db: duckdb.DuckDBPyConnection) -> None:
        db.sql("""
            INSERT INTO variants(
                vid, chrom, pos, ref, alt
            ) VALUES (
                1, 'chr1', 100, 'G', 'C'
            ), (
                2, 'chr2', 200, 'A', 'T'
            );
            
            INSERT INTO vals_info(
                vid, kind, k, v_integer, v_float
            ) VALUES (
                1, 'val', 'af', NULL, 0.3
            ), (
                1, 'val', 'dp', 20, NULL
            ), (
                2, 'val', 'af', NULL, 0.3
            ), (
                2, 'val', 'dp', 20, NULL
            );
            
            INSERT INTO vals_info(
                vid, kind, k, v_varchar
            ) VALUES (
                1, 'info', 'oncokb_oncogenic', 'Oncogenic'
            ), (
                2,'info',  'oncokb_oncogenic', 'not'
            );
        """)

        make_views(db)

        observed = get_somatic_variants_as_df(db)

        expected = pd.DataFrame(
            [
                {
                    "chrom": "chr1",
                    "pos": 100,
                    "ref": "G",
                    "alt": "C",
                    "af": 0.3,
                    "alt_count": None,
                    "am_class": None,
                    "am_pathogenicity": None,
                    "brca1_func_score": None,
                    "civic_description": None,
                    "civic_id": None,
                    "civic_score": None,
                    "cosmic_tier": None,
                    "dbsnp_rs_id": None,
                    "dna_change": None,
                    "dp": 20,
                    "ensembl_feature_id": None,
                    "ensembl_gene_id": None,
                    "exon": None,
                    "gc_content": None,
                    "gnomade_af": None,
                    "gnomadg_af": None,
                    "gt": None,
                    "gtex_gene": None,
                    "gwas_disease": None,
                    "gwas_pmid": None,
                    "hess_driver": False,
                    "hess_signature": None,
                    "hgnc_family": None,
                    "hgnc_name": None,
                    "hugo_symbol": None,
                    "intron": None,
                    "lof": None,
                    "molecular_consequence": None,
                    "nmd": None,
                    "oncogene_high_impact": False,
                    "oncokb_effect": None,
                    "oncokb_hotspot": False,
                    "oncokb_oncogenic": "Oncogenic",
                    "pharmgkb_id": None,
                    "polyphen": None,
                    "protein_change": None,
                    "provean_prediction": None,
                    "ps": None,
                    "ref_count": None,
                    "rescue": True,
                    "revel_score": None,
                    "sift": None,
                    "transcript_likely_lof": None,
                    "tumor_suppressor_high_impact": False,
                    "uniprot_id": None,
                    "variant_info": None,
                    "variant_type": None,
                    "vep_biotype": None,
                    "vep_clin_sig": None,
                    "vep_ensp": None,
                    "vep_existing_variation": None,
                    "vep_hgnc_id": None,
                    "vep_impact": None,
                    "vep_loftool": None,
                    "vep_mane_select": None,
                    "vep_pli_gene_value": None,
                    "vep_somatic": None,
                    "vep_swissprot": None,
                },
            ]
        ).astype(maf_dtypes)

        assert_frame_equal(observed, expected)

    def test_oncokb_hotspot(self, db: duckdb.DuckDBPyConnection) -> None:
        db.sql("""
            INSERT INTO variants(
                vid, chrom, pos, ref, alt
            ) VALUES (
                1, 'chr1', 100, 'G', 'C'
            ), (
                2, 'chr2', 200, 'A', 'T'
            );
            
            INSERT INTO vals_info(
                vid, kind, k, v_integer, v_float
            ) VALUES (
                1, 'val', 'af', NULL, 0.3
            ), (
                1, 'val', 'dp', 20, NULL
            ), (
                2, 'val', 'af', NULL, 0.3
            ), (
                2, 'val', 'dp', 20, NULL
            );
            
            INSERT INTO vals_info(
                vid, kind, k, v_boolean
            ) VALUES (
                1, 'info', 'oncokb_hotspot', TRUE
            );
        """)

        make_views(db)

        observed = get_somatic_variants_as_df(db)

        expected = pd.DataFrame(
            [
                {
                    "chrom": "chr1",
                    "pos": 100,
                    "ref": "G",
                    "alt": "C",
                    "af": 0.3,
                    "alt_count": None,
                    "am_class": None,
                    "am_pathogenicity": None,
                    "brca1_func_score": None,
                    "civic_description": None,
                    "civic_id": None,
                    "civic_score": None,
                    "cosmic_tier": None,
                    "dbsnp_rs_id": None,
                    "dna_change": None,
                    "dp": 20,
                    "ensembl_feature_id": None,
                    "ensembl_gene_id": None,
                    "exon": None,
                    "gc_content": None,
                    "gnomade_af": None,
                    "gnomadg_af": None,
                    "gt": None,
                    "gtex_gene": None,
                    "gwas_disease": None,
                    "gwas_pmid": None,
                    "hess_driver": False,
                    "hess_signature": None,
                    "hgnc_family": None,
                    "hgnc_name": None,
                    "hugo_symbol": None,
                    "intron": None,
                    "lof": None,
                    "molecular_consequence": None,
                    "nmd": None,
                    "oncogene_high_impact": False,
                    "oncokb_effect": None,
                    "oncokb_hotspot": True,
                    "oncokb_oncogenic": None,
                    "pharmgkb_id": None,
                    "polyphen": None,
                    "protein_change": None,
                    "provean_prediction": None,
                    "ps": None,
                    "ref_count": None,
                    "rescue": True,
                    "revel_score": None,
                    "sift": None,
                    "transcript_likely_lof": None,
                    "tumor_suppressor_high_impact": False,
                    "uniprot_id": None,
                    "variant_info": None,
                    "variant_type": None,
                    "vep_biotype": None,
                    "vep_clin_sig": None,
                    "vep_ensp": None,
                    "vep_existing_variation": None,
                    "vep_hgnc_id": None,
                    "vep_impact": None,
                    "vep_loftool": None,
                    "vep_mane_select": None,
                    "vep_pli_gene_value": None,
                    "vep_somatic": None,
                    "vep_swissprot": None,
                },
            ]
        ).astype(maf_dtypes)

        assert_frame_equal(observed, expected)

    def test_cmc_tier(self, db: duckdb.DuckDBPyConnection) -> None:
        db.sql("""
            INSERT INTO variants(
                vid, chrom, pos, ref, alt
            ) VALUES (
                1, 'chr1', 100, 'G', 'C'
            ), (
                2, 'chr2', 200, 'A', 'T'
            );
            
            INSERT INTO vals_info(
                vid, kind, k, v_integer, v_float
            ) VALUES (
                1, 'val', 'af', NULL, 0.3
            ), (
                1, 'val', 'dp', 20, NULL
            ), (
                2, 'val', 'af', NULL, 0.3
            ), (
                2, 'val', 'dp', 20, NULL
            );
            
            INSERT INTO vals_info(
                vid, kind, k, v_integer
            ) VALUES (
                1, 'info', 'cmc_tier', 1 
            ), (
                2, 'info', 'cmc_tier', 3 
            );
        """)

        make_views(db)

        observed = get_somatic_variants_as_df(db)

        expected = pd.DataFrame(
            [
                {
                    "chrom": "chr1",
                    "pos": 100,
                    "ref": "G",
                    "alt": "C",
                    "af": 0.3,
                    "alt_count": None,
                    "am_class": None,
                    "am_pathogenicity": None,
                    "brca1_func_score": None,
                    "civic_description": None,
                    "civic_id": None,
                    "civic_score": None,
                    "cosmic_tier": 1,
                    "dbsnp_rs_id": None,
                    "dna_change": None,
                    "dp": 20,
                    "ensembl_feature_id": None,
                    "ensembl_gene_id": None,
                    "exon": None,
                    "gc_content": None,
                    "gnomade_af": None,
                    "gnomadg_af": None,
                    "gt": None,
                    "gtex_gene": None,
                    "gwas_disease": None,
                    "gwas_pmid": None,
                    "hess_driver": False,
                    "hess_signature": None,
                    "hgnc_family": None,
                    "hgnc_name": None,
                    "hugo_symbol": None,
                    "intron": None,
                    "lof": None,
                    "molecular_consequence": None,
                    "nmd": None,
                    "oncogene_high_impact": False,
                    "oncokb_effect": None,
                    "oncokb_hotspot": False,
                    "oncokb_oncogenic": None,
                    "pharmgkb_id": None,
                    "polyphen": None,
                    "protein_change": None,
                    "provean_prediction": None,
                    "ps": None,
                    "ref_count": None,
                    "rescue": True,
                    "revel_score": None,
                    "sift": None,
                    "transcript_likely_lof": None,
                    "tumor_suppressor_high_impact": False,
                    "uniprot_id": None,
                    "variant_info": None,
                    "variant_type": None,
                    "vep_biotype": None,
                    "vep_clin_sig": None,
                    "vep_ensp": None,
                    "vep_existing_variation": None,
                    "vep_hgnc_id": None,
                    "vep_impact": None,
                    "vep_loftool": None,
                    "vep_mane_select": None,
                    "vep_pli_gene_value": None,
                    "vep_somatic": None,
                    "vep_swissprot": None,
                },
            ]
        ).astype(maf_dtypes)

        assert_frame_equal(observed, expected)

    def test_oc_brca1_func_assay_score(self, db: duckdb.DuckDBPyConnection) -> None:
        db.sql("""
            INSERT INTO variants(
                vid, chrom, pos, ref, alt
            ) VALUES (
                1, 'chr1', 100, 'G', 'C'
            ), (
                2, 'chr2', 200, 'A', 'T'
            );
            
            INSERT INTO vals_info(
                vid, kind, k, v_integer, v_float
            ) VALUES (
                1, 'val', 'af', NULL, 0.3
            ), (
                1, 'val', 'dp', 20, NULL
            ), (
                2, 'val', 'af', NULL, 0.3
            ), (
                2, 'val', 'dp', 20, NULL
            );
            
            INSERT INTO vals_info(
                vid, kind, k, v_float
            ) VALUES (
                1, 'info', 'oc_brca1_func_assay_score', -1.5
            ), (
                2, 'info', 'oc_brca1_func_assay_score', -0.9
            );
        """)

        make_views(db, max_brca1_func_assay_score=-1.328)

        observed = get_somatic_variants_as_df(db)

        expected = pd.DataFrame(
            [
                {
                    "chrom": "chr1",
                    "pos": 100,
                    "ref": "G",
                    "alt": "C",
                    "af": 0.3,
                    "alt_count": None,
                    "am_class": None,
                    "am_pathogenicity": None,
                    "brca1_func_score": -1.5,
                    "civic_description": None,
                    "civic_id": None,
                    "civic_score": None,
                    "cosmic_tier": None,
                    "dbsnp_rs_id": None,
                    "dna_change": None,
                    "dp": 20,
                    "ensembl_feature_id": None,
                    "ensembl_gene_id": None,
                    "exon": None,
                    "gc_content": None,
                    "gnomade_af": None,
                    "gnomadg_af": None,
                    "gt": None,
                    "gtex_gene": None,
                    "gwas_disease": None,
                    "gwas_pmid": None,
                    "hess_driver": False,
                    "hess_signature": None,
                    "hgnc_family": None,
                    "hgnc_name": None,
                    "hugo_symbol": None,
                    "intron": None,
                    "lof": None,
                    "molecular_consequence": None,
                    "nmd": None,
                    "oncogene_high_impact": False,
                    "oncokb_effect": None,
                    "oncokb_hotspot": False,
                    "oncokb_oncogenic": None,
                    "pharmgkb_id": None,
                    "polyphen": None,
                    "protein_change": None,
                    "provean_prediction": None,
                    "ps": None,
                    "ref_count": None,
                    "rescue": True,
                    "revel_score": None,
                    "sift": None,
                    "transcript_likely_lof": None,
                    "tumor_suppressor_high_impact": False,
                    "uniprot_id": None,
                    "variant_info": None,
                    "variant_type": None,
                    "vep_biotype": None,
                    "vep_clin_sig": None,
                    "vep_ensp": None,
                    "vep_existing_variation": None,
                    "vep_hgnc_id": None,
                    "vep_impact": None,
                    "vep_loftool": None,
                    "vep_mane_select": None,
                    "vep_pli_gene_value": None,
                    "vep_somatic": None,
                    "vep_swissprot": None,
                },
            ]
        ).astype(maf_dtypes)

        assert_frame_equal(observed, expected)

    def test_oncogene(self, db: duckdb.DuckDBPyConnection) -> None:
        db.sql("""
            INSERT INTO variants(
                vid, chrom, pos, ref, alt
            ) VALUES (
                1, 'chr1', 100, 'G', 'C'
            ), (
                2, 'chr2', 200, 'A', 'T'
            ), (
                3, 'chr3', 300, 'A', 'T'
            );
            
            INSERT INTO vals_info(
                vid, kind, k, v_integer, v_float
            ) VALUES (
                1, 'val', 'af', NULL, 0.3
            ), (
                1, 'val', 'dp', 20, NULL
            ), (
                2, 'val', 'af', NULL, 0.3
            ), (
                2, 'val', 'dp', 20, NULL
            ), (
                3, 'val', 'af', NULL, 0.3
            ), (
                3, 'val', 'dp', 20, NULL
            );
            
            INSERT INTO vals_info(
                vid, kind, k, v_boolean, v_json_arr
            ) VALUES (
                1, 'info', 'oncogene', TRUE, NULL 
            ), (
                1, 'info', 'csq', NULL, ['{"impact": "HIGH"}']
            ), (
                2, 'info', 'oncogene', TRUE, NULL 
            ), (
                2, 'info', 'csq', NULL, ['{"impact": "LOW"}']
            ), (
                3, 'info', 'oncogene', TRUE, NULL 
            );
        """)

        make_views(db)

        observed = get_somatic_variants_as_df(db)

        expected = pd.DataFrame(
            [
                {
                    "chrom": "chr1",
                    "pos": 100,
                    "ref": "G",
                    "alt": "C",
                    "af": 0.3,
                    "alt_count": None,
                    "am_class": None,
                    "am_pathogenicity": None,
                    "brca1_func_score": None,
                    "civic_description": None,
                    "civic_id": None,
                    "civic_score": None,
                    "cosmic_tier": None,
                    "dbsnp_rs_id": None,
                    "dna_change": None,
                    "dp": 20,
                    "ensembl_feature_id": None,
                    "ensembl_gene_id": None,
                    "exon": None,
                    "gc_content": None,
                    "gnomade_af": None,
                    "gnomadg_af": None,
                    "gt": None,
                    "gtex_gene": None,
                    "gwas_disease": None,
                    "gwas_pmid": None,
                    "hess_driver": False,
                    "hess_signature": None,
                    "hgnc_family": None,
                    "hgnc_name": None,
                    "hugo_symbol": None,
                    "intron": None,
                    "lof": None,
                    "molecular_consequence": None,
                    "nmd": None,
                    "oncogene_high_impact": True,
                    "oncokb_effect": None,
                    "oncokb_hotspot": False,
                    "oncokb_oncogenic": None,
                    "pharmgkb_id": None,
                    "polyphen": None,
                    "protein_change": None,
                    "provean_prediction": None,
                    "ps": None,
                    "ref_count": None,
                    "rescue": True,
                    "revel_score": None,
                    "sift": None,
                    "transcript_likely_lof": None,
                    "tumor_suppressor_high_impact": False,
                    "uniprot_id": None,
                    "variant_info": None,
                    "variant_type": None,
                    "vep_biotype": None,
                    "vep_clin_sig": None,
                    "vep_ensp": None,
                    "vep_existing_variation": None,
                    "vep_hgnc_id": None,
                    "vep_impact": "HIGH",
                    "vep_loftool": None,
                    "vep_mane_select": None,
                    "vep_pli_gene_value": None,
                    "vep_somatic": None,
                    "vep_swissprot": None,
                },
            ]
        ).astype(maf_dtypes)

        assert_frame_equal(observed, expected)

    def test_tsg(self, db: duckdb.DuckDBPyConnection) -> None:
        db.sql("""
            INSERT INTO variants(
                vid, chrom, pos, ref, alt
            ) VALUES (
                1, 'chr1', 100, 'G', 'C'
            ), (
                2, 'chr2', 200, 'A', 'T'
            ), (
                3, 'chr3', 300, 'A', 'T'
            );
            
            INSERT INTO vals_info(
                vid, kind, k, v_integer, v_float
            ) VALUES (
                1, 'val', 'af', NULL, 0.3
            ), (
                1, 'val', 'dp', 20, NULL
            ), (
                2, 'val', 'af', NULL, 0.3
            ), (
                2, 'val', 'dp', 20, NULL
            ), (
                3, 'val', 'af', NULL, 0.3
            ), (
                3, 'val', 'dp', 20, NULL
            );
            
            INSERT INTO vals_info(
                vid, kind, k, v_boolean, v_json_arr
            ) VALUES (
                1, 'info', 'tsg', TRUE, NULL 
            ), (
                1, 'info', 'csq', NULL, ['{"impact": "HIGH"}']
            ), (
                2, 'info', 'tsg', TRUE, NULL 
            ), (
                2, 'info', 'csq', NULL, ['{"impact": "LOW"}']
            ), (
                3, 'info', 'tsg', TRUE, NULL 
            );
        """)

        make_views(db)

        observed = get_somatic_variants_as_df(db)

        expected = pd.DataFrame(
            [
                {
                    "chrom": "chr1",
                    "pos": 100,
                    "ref": "G",
                    "alt": "C",
                    "af": 0.3,
                    "alt_count": None,
                    "am_class": None,
                    "am_pathogenicity": None,
                    "brca1_func_score": None,
                    "civic_description": None,
                    "civic_id": None,
                    "civic_score": None,
                    "cosmic_tier": None,
                    "dbsnp_rs_id": None,
                    "dna_change": None,
                    "dp": 20,
                    "ensembl_feature_id": None,
                    "ensembl_gene_id": None,
                    "exon": None,
                    "gc_content": None,
                    "gnomade_af": None,
                    "gnomadg_af": None,
                    "gt": None,
                    "gtex_gene": None,
                    "gwas_disease": None,
                    "gwas_pmid": None,
                    "hess_driver": False,
                    "hess_signature": None,
                    "hgnc_family": None,
                    "hgnc_name": None,
                    "hugo_symbol": None,
                    "intron": None,
                    "lof": None,
                    "molecular_consequence": None,
                    "nmd": None,
                    "oncogene_high_impact": False,
                    "oncokb_effect": None,
                    "oncokb_hotspot": False,
                    "oncokb_oncogenic": None,
                    "pharmgkb_id": None,
                    "polyphen": None,
                    "protein_change": None,
                    "provean_prediction": None,
                    "ps": None,
                    "ref_count": None,
                    "rescue": True,
                    "revel_score": None,
                    "sift": None,
                    "transcript_likely_lof": None,
                    "tumor_suppressor_high_impact": True,
                    "uniprot_id": None,
                    "variant_info": None,
                    "variant_type": None,
                    "vep_biotype": None,
                    "vep_clin_sig": None,
                    "vep_ensp": None,
                    "vep_existing_variation": None,
                    "vep_hgnc_id": None,
                    "vep_impact": "HIGH",
                    "vep_loftool": None,
                    "vep_mane_select": None,
                    "vep_pli_gene_value": None,
                    "vep_somatic": None,
                    "vep_swissprot": None,
                },
            ]
        ).astype(maf_dtypes)

        assert_frame_equal(observed, expected)

    def test_hess_driver(self, db: duckdb.DuckDBPyConnection) -> None:
        db.sql("""
            INSERT INTO variants(
                vid, chrom, pos, ref, alt
            ) VALUES (
                1, 'chr1', 100, 'G', 'C'
            ), (
                2, 'chr2', 200, 'A', 'T'
            );
            
            INSERT INTO vals_info(
                vid, kind, k, v_integer, v_float
            ) VALUES (
                1, 'val', 'af', NULL, 0.3
            ), (
                1, 'val', 'dp', 20, NULL
            ), (
                2, 'val', 'af', NULL, 0.3
            ), (
                2, 'val', 'dp', 20, NULL
            );
            
            INSERT INTO vals_info(
                vid, kind, k, v_varchar
            ) VALUES (
                1, 'info', 'hess', 'sig'
            );
        """)

        make_views(db)

        observed = get_somatic_variants_as_df(db)

        expected = pd.DataFrame(
            [
                {
                    "chrom": "chr1",
                    "pos": 100,
                    "ref": "G",
                    "alt": "C",
                    "af": 0.3,
                    "alt_count": None,
                    "am_class": None,
                    "am_pathogenicity": None,
                    "brca1_func_score": None,
                    "civic_description": None,
                    "civic_id": None,
                    "civic_score": None,
                    "cosmic_tier": None,
                    "dbsnp_rs_id": None,
                    "dna_change": None,
                    "dp": 20,
                    "ensembl_feature_id": None,
                    "ensembl_gene_id": None,
                    "exon": None,
                    "gc_content": None,
                    "gnomade_af": None,
                    "gnomadg_af": None,
                    "gt": None,
                    "gtex_gene": None,
                    "gwas_disease": None,
                    "gwas_pmid": None,
                    "hess_driver": True,
                    "hess_signature": "sig",
                    "hgnc_family": None,
                    "hgnc_name": None,
                    "hugo_symbol": None,
                    "intron": None,
                    "lof": None,
                    "molecular_consequence": None,
                    "nmd": None,
                    "oncogene_high_impact": False,
                    "oncokb_effect": None,
                    "oncokb_hotspot": False,
                    "oncokb_oncogenic": None,
                    "pharmgkb_id": None,
                    "polyphen": None,
                    "protein_change": None,
                    "provean_prediction": None,
                    "ps": None,
                    "ref_count": None,
                    "rescue": True,
                    "revel_score": None,
                    "sift": None,
                    "transcript_likely_lof": None,
                    "tumor_suppressor_high_impact": False,
                    "uniprot_id": None,
                    "variant_info": None,
                    "variant_type": None,
                    "vep_biotype": None,
                    "vep_clin_sig": None,
                    "vep_ensp": None,
                    "vep_existing_variation": None,
                    "vep_hgnc_id": None,
                    "vep_impact": None,
                    "vep_loftool": None,
                    "vep_mane_select": None,
                    "vep_pli_gene_value": None,
                    "vep_somatic": None,
                    "vep_swissprot": None,
                },
            ]
        ).astype(maf_dtypes)

        assert_frame_equal(observed, expected)

    def test_tert(self, db: duckdb.DuckDBPyConnection) -> None:
        db.sql("""
            INSERT INTO variants(
                vid, chrom, pos, ref, alt
            ) VALUES (
                1, 'chr5', 1295059, 'G', 'C'
            ), (
                2, 'chr2', 200, 'A', 'T'
            );
            
            INSERT INTO vals_info(
                vid, kind, k, v_integer, v_float
            ) VALUES (
                1, 'val', 'af', NULL, 0.3
            ), (
                1, 'val', 'dp', 20, NULL
            ), (
                2, 'val', 'af', NULL, 0.3
            ), (
                2, 'val', 'dp', 20, NULL
            );
            
            INSERT INTO vals_info(
                vid, kind, k, v_json_arr
            ) VALUES (
                1, 'info', 'csq', ['{"symbol": "TERT"}']
            );
        """)

        make_views(db)

        observed = get_somatic_variants_as_df(db)

        expected = pd.DataFrame(
            [
                {
                    "chrom": "chr5",
                    "pos": 1295059,
                    "ref": "G",
                    "alt": "C",
                    "af": 0.3,
                    "alt_count": None,
                    "am_class": None,
                    "am_pathogenicity": None,
                    "brca1_func_score": None,
                    "civic_description": None,
                    "civic_id": None,
                    "civic_score": None,
                    "cosmic_tier": None,
                    "dbsnp_rs_id": None,
                    "dna_change": None,
                    "dp": 20,
                    "ensembl_feature_id": None,
                    "ensembl_gene_id": None,
                    "exon": None,
                    "gc_content": None,
                    "gnomade_af": None,
                    "gnomadg_af": None,
                    "gt": None,
                    "gtex_gene": None,
                    "gwas_disease": None,
                    "gwas_pmid": None,
                    "hess_driver": False,
                    "hess_signature": None,
                    "hgnc_family": None,
                    "hgnc_name": None,
                    "hugo_symbol": "TERT",
                    "intron": None,
                    "lof": None,
                    "molecular_consequence": None,
                    "nmd": None,
                    "oncogene_high_impact": False,
                    "oncokb_effect": None,
                    "oncokb_hotspot": False,
                    "oncokb_oncogenic": None,
                    "pharmgkb_id": None,
                    "polyphen": None,
                    "protein_change": None,
                    "provean_prediction": None,
                    "ps": None,
                    "ref_count": None,
                    "rescue": True,
                    "revel_score": None,
                    "sift": None,
                    "transcript_likely_lof": None,
                    "tumor_suppressor_high_impact": False,
                    "uniprot_id": None,
                    "variant_info": None,
                    "variant_type": None,
                    "vep_biotype": None,
                    "vep_clin_sig": None,
                    "vep_ensp": None,
                    "vep_existing_variation": None,
                    "vep_hgnc_id": None,
                    "vep_impact": None,
                    "vep_loftool": None,
                    "vep_mane_select": None,
                    "vep_pli_gene_value": None,
                    "vep_somatic": None,
                    "vep_swissprot": None,
                },
            ]
        ).astype(maf_dtypes)

        assert_frame_equal(observed, expected)

    def test_met(self, db: duckdb.DuckDBPyConnection) -> None:
        db.sql("""
            INSERT INTO variants(
                vid, chrom, pos, ref, alt
            ) VALUES (
                1, 'chr7', 116771829, 'G', 'C'
            ), (
                2, 'chr2', 200, 'A', 'T'
            );
            
            INSERT INTO vals_info(
                vid, kind, k, v_integer, v_float
            ) VALUES (
                1, 'val', 'af', NULL, 0.3
            ), (
                1, 'val', 'dp', 20, NULL
            ), (
                2, 'val', 'af', NULL, 0.3
            ), (
                2, 'val', 'dp', 20, NULL
            );
            
            INSERT INTO vals_info(
                vid, kind, k, v_json_arr
            ) VALUES (
                1, 'info', 'csq', ['{"symbol": "MET"}']
            );
        """)

        make_views(db)

        observed = get_somatic_variants_as_df(db)

        expected = pd.DataFrame(
            [
                {
                    "chrom": "chr7",
                    "pos": 116771829,
                    "ref": "G",
                    "alt": "C",
                    "af": 0.3,
                    "alt_count": None,
                    "am_class": None,
                    "am_pathogenicity": None,
                    "brca1_func_score": None,
                    "civic_description": None,
                    "civic_id": None,
                    "civic_score": None,
                    "cosmic_tier": None,
                    "dbsnp_rs_id": None,
                    "dna_change": None,
                    "dp": 20,
                    "ensembl_feature_id": None,
                    "ensembl_gene_id": None,
                    "exon": None,
                    "gc_content": None,
                    "gnomade_af": None,
                    "gnomadg_af": None,
                    "gt": None,
                    "gtex_gene": None,
                    "gwas_disease": None,
                    "gwas_pmid": None,
                    "hess_driver": False,
                    "hess_signature": None,
                    "hgnc_family": None,
                    "hgnc_name": None,
                    "hugo_symbol": "MET",
                    "intron": None,
                    "lof": None,
                    "molecular_consequence": None,
                    "nmd": None,
                    "oncogene_high_impact": False,
                    "oncokb_effect": None,
                    "oncokb_hotspot": False,
                    "oncokb_oncogenic": None,
                    "pharmgkb_id": None,
                    "polyphen": None,
                    "protein_change": None,
                    "provean_prediction": None,
                    "ps": None,
                    "ref_count": None,
                    "rescue": True,
                    "revel_score": None,
                    "sift": None,
                    "transcript_likely_lof": None,
                    "tumor_suppressor_high_impact": False,
                    "uniprot_id": None,
                    "variant_info": None,
                    "variant_type": None,
                    "vep_biotype": None,
                    "vep_clin_sig": None,
                    "vep_ensp": None,
                    "vep_existing_variation": None,
                    "vep_hgnc_id": None,
                    "vep_impact": None,
                    "vep_loftool": None,
                    "vep_mane_select": None,
                    "vep_pli_gene_value": None,
                    "vep_somatic": None,
                    "vep_swissprot": None,
                },
            ]
        ).astype(maf_dtypes)

        assert_frame_equal(observed, expected)


@pytest.mark.usefixtures("db_setup")
class TestFilteredVids:
    def test_splice_event(self, db: duckdb.DuckDBPyConnection) -> None:
        db.sql("""
            INSERT INTO variants(
                vid, chrom, pos, ref, alt
            ) VALUES (
                1, 'chr1', 100, 'G', 'C'
            ), (
                2, 'chr2', 200, 'A', 'T'
            ), (
                3, 'chr3', 300, 'A', 'T'
            );
            
            INSERT INTO vals_info(
                vid, kind, k, v_integer, v_float
            ) VALUES (
                1, 'val', 'af', NULL, 0.3
            ), (
                1, 'val', 'dp', 20, NULL
            ), (
                2, 'val', 'af', NULL, 0.3
            ), (
                2, 'val', 'dp', 20, NULL
            ), (
                3, 'val', 'af', NULL, 0.3
            ), (
                3, 'val', 'dp', 20, NULL
            );
            
            INSERT INTO vals_info(
                vid, kind, k, v_json_arr
            ) VALUES (
                1, 'info', 'csq', ['{"impact": "HIGH", "consequence": "foo&splice_region"}']
            ), (
                2, 'info', 'csq', ['{"impact": "HIGH", "consequence": "other"}']
            ), (
                3, 'info', 'csq', ['{"impact": "LOW", "consequence": "splice"}']
            );
        """)

        make_views(db)

        observed = get_somatic_variants_as_df(db)

        expected = pd.DataFrame(
            [
                {
                    "chrom": "chr1",
                    "pos": 100,
                    "ref": "G",
                    "alt": "C",
                    "af": 0.3,
                    "alt_count": None,
                    "am_class": None,
                    "am_pathogenicity": None,
                    "brca1_func_score": None,
                    "civic_description": None,
                    "civic_id": None,
                    "civic_score": None,
                    "cosmic_tier": None,
                    "dbsnp_rs_id": None,
                    "dna_change": None,
                    "dp": 20,
                    "ensembl_feature_id": None,
                    "ensembl_gene_id": None,
                    "exon": None,
                    "gc_content": None,
                    "gnomade_af": None,
                    "gnomadg_af": None,
                    "gt": None,
                    "gtex_gene": None,
                    "gwas_disease": None,
                    "gwas_pmid": None,
                    "hess_driver": False,
                    "hess_signature": None,
                    "hgnc_family": None,
                    "hgnc_name": None,
                    "hugo_symbol": None,
                    "intron": None,
                    "lof": None,
                    "molecular_consequence": None,
                    "nmd": None,
                    "oncogene_high_impact": False,
                    "oncokb_effect": None,
                    "oncokb_hotspot": False,
                    "oncokb_oncogenic": None,
                    "pharmgkb_id": None,
                    "polyphen": None,
                    "protein_change": None,
                    "provean_prediction": None,
                    "ps": None,
                    "ref_count": None,
                    "rescue": False,
                    "revel_score": None,
                    "sift": None,
                    "transcript_likely_lof": None,
                    "tumor_suppressor_high_impact": False,
                    "uniprot_id": None,
                    "variant_info": "foo&splice_region",
                    "variant_type": None,
                    "vep_biotype": None,
                    "vep_clin_sig": None,
                    "vep_ensp": None,
                    "vep_existing_variation": None,
                    "vep_hgnc_id": None,
                    "vep_impact": "HIGH",
                    "vep_loftool": None,
                    "vep_mane_select": None,
                    "vep_pli_gene_value": None,
                    "vep_somatic": None,
                    "vep_swissprot": None,
                },
            ]
        ).astype(maf_dtypes)

        assert_frame_equal(observed, expected)

    def test_protein_change(self, db: duckdb.DuckDBPyConnection) -> None:
        db.sql("""
            INSERT INTO variants(
                vid, chrom, pos, ref, alt
            ) VALUES (
                1, 'chr1', 100, 'G', 'C'
            ), (
                2, 'chr2', 200, 'A', 'T'
            ), (
                3, 'chr3', 300, 'A', 'T'
            );
            
            INSERT INTO vals_info(
                vid, kind, k, v_integer, v_float
            ) VALUES (
                1, 'val', 'af', NULL, 0.3
            ), (
                1, 'val', 'dp', 20, NULL
            ), (
                2, 'val', 'af', NULL, 0.3
            ), (
                2, 'val', 'dp', 20, NULL
            ), (
                3, 'val', 'af', NULL, 0.3
            ), (
                3, 'val', 'dp', 20, NULL
            );
            
            INSERT INTO vals_info(
                vid, kind, k, v_json_arr
            ) VALUES (
                1, 'info', 'csq', ['{"impact": "HIGH", "hgvsp": "ENSP15.3:p.Glu32del"}']
            ), (
                2, 'info', 'csq', ['{"impact": "HIGH", "hgvsp": "ENSP30.3:p.Asp637="}']
            ), (
                3, 'info', 'csq', ['{"impact": "LOW", "hgvsp": null}']
            );
        """)

        make_views(db)

        observed = get_somatic_variants_as_df(db)

        expected = pd.DataFrame(
            [
                {
                    "chrom": "chr1",
                    "pos": 100,
                    "ref": "G",
                    "alt": "C",
                    "af": 0.3,
                    "alt_count": None,
                    "am_class": None,
                    "am_pathogenicity": None,
                    "brca1_func_score": None,
                    "civic_description": None,
                    "civic_id": None,
                    "civic_score": None,
                    "cosmic_tier": None,
                    "dbsnp_rs_id": None,
                    "dna_change": None,
                    "dp": 20,
                    "ensembl_feature_id": None,
                    "ensembl_gene_id": None,
                    "exon": None,
                    "gc_content": None,
                    "gnomade_af": None,
                    "gnomadg_af": None,
                    "gt": None,
                    "gtex_gene": None,
                    "gwas_disease": None,
                    "gwas_pmid": None,
                    "hess_driver": False,
                    "hess_signature": None,
                    "hgnc_family": None,
                    "hgnc_name": None,
                    "hugo_symbol": None,
                    "intron": None,
                    "lof": None,
                    "molecular_consequence": None,
                    "nmd": None,
                    "oncogene_high_impact": False,
                    "oncokb_effect": None,
                    "oncokb_hotspot": False,
                    "oncokb_oncogenic": None,
                    "pharmgkb_id": None,
                    "polyphen": None,
                    "protein_change": "ENSP15.3:p.Glu32del",
                    "provean_prediction": None,
                    "ps": None,
                    "ref_count": None,
                    "rescue": False,
                    "revel_score": None,
                    "sift": None,
                    "transcript_likely_lof": None,
                    "tumor_suppressor_high_impact": False,
                    "uniprot_id": None,
                    "variant_info": None,
                    "variant_type": None,
                    "vep_biotype": None,
                    "vep_clin_sig": None,
                    "vep_ensp": None,
                    "vep_existing_variation": None,
                    "vep_hgnc_id": None,
                    "vep_impact": "HIGH",
                    "vep_loftool": None,
                    "vep_mane_select": None,
                    "vep_pli_gene_value": None,
                    "vep_somatic": None,
                    "vep_swissprot": None,
                },
            ]
        ).astype(maf_dtypes)

        assert_frame_equal(observed, expected)

    def test_clustered_event(self, db: duckdb.DuckDBPyConnection) -> None:
        db.sql("""
            INSERT INTO variants(
                vid, chrom, pos, ref, alt, filters
            ) VALUES (
                1, 'chr1', 100, 'G', 'C', ['germline']
            ), (
                2, 'chr2', 200, 'A', 'T', ['clustered_event']
            );
            
            INSERT INTO vals_info(
                vid, kind, k, v_integer, v_float
            ) VALUES (
                1, 'val', 'af', NULL, 0.3
            ), (
                1, 'val', 'dp', 20, NULL
            ), (
                2, 'val', 'af', NULL, 0.3
            ), (
                2, 'val', 'dp', 20, NULL
            );
            
            INSERT INTO vals_info(
                vid, kind, k, v_json_arr
            ) VALUES (
                1, 'info', 'csq', ['{"impact": "HIGH", "consequence": "splice", "gnom_ade_af": 0.000001, "gnom_adg_af": 0.000002}']
            ), (
                2, 'info', 'csq', ['{"impact": "HIGH", "consequence": "splice", "gnom_ade_af": 0.000001, "gnom_adg_af": 0.000002}']
            );
        """)

        make_views(db)

        observed = get_somatic_variants_as_df(db)

        expected = pd.DataFrame(
            [
                {
                    "chrom": "chr1",
                    "pos": 100,
                    "ref": "G",
                    "alt": "C",
                    "af": 0.3,
                    "alt_count": None,
                    "am_class": None,
                    "am_pathogenicity": None,
                    "brca1_func_score": None,
                    "civic_description": None,
                    "civic_id": None,
                    "civic_score": None,
                    "cosmic_tier": None,
                    "dbsnp_rs_id": None,
                    "dna_change": None,
                    "dp": 20,
                    "ensembl_feature_id": None,
                    "ensembl_gene_id": None,
                    "exon": None,
                    "gc_content": None,
                    "gnomade_af": 0.000001,
                    "gnomadg_af": 0.000002,
                    "gt": None,
                    "gtex_gene": None,
                    "gwas_disease": None,
                    "gwas_pmid": None,
                    "hess_driver": False,
                    "hess_signature": None,
                    "hgnc_family": None,
                    "hgnc_name": None,
                    "hugo_symbol": None,
                    "intron": None,
                    "lof": None,
                    "molecular_consequence": None,
                    "nmd": None,
                    "oncogene_high_impact": False,
                    "oncokb_effect": None,
                    "oncokb_hotspot": False,
                    "oncokb_oncogenic": None,
                    "pharmgkb_id": None,
                    "polyphen": None,
                    "protein_change": None,
                    "provean_prediction": None,
                    "ps": None,
                    "ref_count": None,
                    "rescue": False,
                    "revel_score": None,
                    "sift": None,
                    "transcript_likely_lof": None,
                    "tumor_suppressor_high_impact": False,
                    "uniprot_id": None,
                    "variant_info": "splice",
                    "variant_type": None,
                    "vep_biotype": None,
                    "vep_clin_sig": None,
                    "vep_ensp": None,
                    "vep_existing_variation": None,
                    "vep_hgnc_id": None,
                    "vep_impact": "HIGH",
                    "vep_loftool": None,
                    "vep_mane_select": None,
                    "vep_pli_gene_value": None,
                    "vep_somatic": None,
                    "vep_swissprot": None,
                },
            ]
        ).astype(maf_dtypes)

        assert_frame_equal(observed, expected)

    def test_segdup_rm(self, db: duckdb.DuckDBPyConnection) -> None:
        db.sql("""
            INSERT INTO variants(
                vid, chrom, pos, ref, alt
            ) VALUES (
                1, 'chr1', 100, 'G', 'C'
            ), (
                2, 'chr2', 200, 'A', 'T'
            ), (
                3, 'chr3', 300, 'A', 'T'
            );
            
            INSERT INTO vals_info(
                vid, kind, k, v_integer, v_float
            ) VALUES (
                1, 'val', 'af', NULL, 0.3
            ), (
                1, 'val', 'dp', 20, NULL
            ), (
                2, 'val', 'af', NULL, 0.3
            ), (
                2, 'val', 'dp', 20, NULL
            ), (
                3, 'val', 'af', NULL, 0.3
            ), (
                3, 'val', 'dp', 20, NULL
            );
            
            INSERT INTO vals_info(
                vid, kind, k, v_json_arr
            ) VALUES (
                1, 'info', 'csq', ['{"impact": "HIGH", "consequence": "splice"}']
            ), (
                2, 'info', 'csq', ['{"impact": "HIGH", "consequence": "splice"}']
            ), (
                3, 'info', 'csq', ['{"impact": "HIGH", "consequence": "splice"}']
            );
            
            --undo some of the previous positive selections
            INSERT INTO vals_info(
                vid, kind, k, v_boolean
            ) VALUES (
                2, 'info', 'segdup', TRUE
            ), (
                3, 'info', 'rm', TRUE
            );
        """)

        make_views(db)

        observed = get_somatic_variants_as_df(db)

        expected = pd.DataFrame(
            [
                {
                    "chrom": "chr1",
                    "pos": 100,
                    "ref": "G",
                    "alt": "C",
                    "af": 0.3,
                    "alt_count": None,
                    "am_class": None,
                    "am_pathogenicity": None,
                    "brca1_func_score": None,
                    "civic_description": None,
                    "civic_id": None,
                    "civic_score": None,
                    "cosmic_tier": None,
                    "dbsnp_rs_id": None,
                    "dna_change": None,
                    "dp": 20,
                    "ensembl_feature_id": None,
                    "ensembl_gene_id": None,
                    "exon": None,
                    "gc_content": None,
                    "gnomade_af": None,
                    "gnomadg_af": None,
                    "gt": None,
                    "gtex_gene": None,
                    "gwas_disease": None,
                    "gwas_pmid": None,
                    "hess_driver": False,
                    "hess_signature": None,
                    "hgnc_family": None,
                    "hgnc_name": None,
                    "hugo_symbol": None,
                    "intron": None,
                    "lof": None,
                    "molecular_consequence": None,
                    "nmd": None,
                    "oncogene_high_impact": False,
                    "oncokb_effect": None,
                    "oncokb_hotspot": False,
                    "oncokb_oncogenic": None,
                    "pharmgkb_id": None,
                    "polyphen": None,
                    "protein_change": None,
                    "provean_prediction": None,
                    "ps": None,
                    "ref_count": None,
                    "rescue": False,
                    "revel_score": None,
                    "sift": None,
                    "transcript_likely_lof": None,
                    "tumor_suppressor_high_impact": False,
                    "uniprot_id": None,
                    "variant_info": "splice",
                    "variant_type": None,
                    "vep_biotype": None,
                    "vep_clin_sig": None,
                    "vep_ensp": None,
                    "vep_existing_variation": None,
                    "vep_hgnc_id": None,
                    "vep_impact": "HIGH",
                    "vep_loftool": None,
                    "vep_mane_select": None,
                    "vep_pli_gene_value": None,
                    "vep_somatic": None,
                    "vep_swissprot": None,
                },
            ]
        ).astype(maf_dtypes)

        assert_frame_equal(observed, expected)

    def test_gnomad(self, db: duckdb.DuckDBPyConnection) -> None:
        db.sql("""
            INSERT INTO variants(
                vid, chrom, pos, ref, alt
            ) VALUES (
                1, 'chr1', 100, 'G', 'C'
            ), (
                2, 'chr2', 200, 'A', 'T'
            ), (
                3, 'chr3', 300, 'A', 'T'
            ), (
                4, 'chr4', 400, 'A', 'T'
            );
            
            INSERT INTO vals_info(
                vid, kind, k, v_integer, v_float
            ) VALUES (
                1, 'val', 'af', NULL, 0.3
            ), (
                1, 'val', 'dp', 20, NULL
            ), (
                2, 'val', 'af', NULL, 0.3
            ), (
                2, 'val', 'dp', 20, NULL
            ), (
                3, 'val', 'af', NULL, 0.3
            ), (
                3, 'val', 'dp', 20, NULL
            ), (
                4, 'val', 'af', NULL, 0.3
            ), (
                4, 'val', 'dp', 20, NULL
            );
            
            INSERT INTO vals_info(
                vid, kind, k, v_json_arr
            ) VALUES (
                1, 'info', 'csq', ['{"impact": "HIGH", "consequence": "splice", "gnom_ade_af": 0.000001, "gnom_adg_af": 0.000002}']
            ), (
                2, 'info', 'csq', ['{"impact": "HIGH", "consequence": "splice", "gnom_ade_af": 0.000001, "gnom_adg_af": 0.1}']
            ), (
                3, 'info', 'csq', ['{"impact": "HIGH", "consequence": "splice", "gnom_ade_af": null, "gnom_adg_af": 0.1}']
            ), (
                4, 'info', 'csq', ['{"impact": "HIGH", "consequence": "splice", "gnom_ade_af": 0.000001, "gnom_adg_af": 0.000002}']
            );
            
            INSERT INTO vals_info(
                vid, kind, k, v_boolean
            ) VALUES (
                4, 'info', 'pon', TRUE
            )
        """)

        make_views(db)

        observed = get_somatic_variants_as_df(db)

        expected = pd.DataFrame(
            [
                {
                    "chrom": "chr1",
                    "pos": 100,
                    "ref": "G",
                    "alt": "C",
                    "af": 0.3,
                    "alt_count": None,
                    "am_class": None,
                    "am_pathogenicity": None,
                    "brca1_func_score": None,
                    "civic_description": None,
                    "civic_id": None,
                    "civic_score": None,
                    "cosmic_tier": None,
                    "dbsnp_rs_id": None,
                    "dna_change": None,
                    "dp": 20,
                    "ensembl_feature_id": None,
                    "ensembl_gene_id": None,
                    "exon": None,
                    "gc_content": None,
                    "gnomade_af": 0.000001,
                    "gnomadg_af": 0.000002,
                    "gt": None,
                    "gtex_gene": None,
                    "gwas_disease": None,
                    "gwas_pmid": None,
                    "hess_driver": False,
                    "hess_signature": None,
                    "hgnc_family": None,
                    "hgnc_name": None,
                    "hugo_symbol": None,
                    "intron": None,
                    "lof": None,
                    "molecular_consequence": None,
                    "nmd": None,
                    "oncogene_high_impact": False,
                    "oncokb_effect": None,
                    "oncokb_hotspot": False,
                    "oncokb_oncogenic": None,
                    "pharmgkb_id": None,
                    "polyphen": None,
                    "protein_change": None,
                    "provean_prediction": None,
                    "ps": None,
                    "ref_count": None,
                    "rescue": False,
                    "revel_score": None,
                    "sift": None,
                    "transcript_likely_lof": None,
                    "tumor_suppressor_high_impact": False,
                    "uniprot_id": None,
                    "variant_info": "splice",
                    "variant_type": None,
                    "vep_biotype": None,
                    "vep_clin_sig": None,
                    "vep_ensp": None,
                    "vep_existing_variation": None,
                    "vep_hgnc_id": None,
                    "vep_impact": "HIGH",
                    "vep_loftool": None,
                    "vep_mane_select": None,
                    "vep_pli_gene_value": None,
                    "vep_somatic": None,
                    "vep_swissprot": None,
                },
            ]
        ).astype(maf_dtypes)

        assert_frame_equal(observed, expected)


@pytest.mark.usefixtures("db_setup")
class TestColumns:
    def test_basic(self, db: duckdb.DuckDBPyConnection) -> None:
        db.sql("""
            INSERT INTO variants(
                vid, chrom, pos, ref, alt
            ) VALUES (
                1, 'chr1', 100, 'G', 'C'
            );
            
            INSERT INTO vals_info(
                vid, kind, k, v_varchar, v_integer, v_float, v_integer_arr
            ) VALUES (
                1, 'val', 'af', NULL, NULL, 0.3, NULL
            ), (
                1, 'val', 'dp', NULL, 20, NULL, NULL
            );
            
            --any rescue state will work for this test
            INSERT INTO vals_info(
                vid, kind, k, v_varchar
            ) VALUES (
                1, 'info', 'oncokb_muteff', 'Loss-of-function'
            );
        """)

        make_views(db)

        observed = get_somatic_variants_as_df(db)

        expected = pd.DataFrame(
            [
                {
                    "chrom": "chr1",
                    "pos": 100,
                    "ref": "G",
                    "alt": "C",
                    "af": 0.3,
                    "alt_count": None,
                    "am_class": None,
                    "am_pathogenicity": None,
                    "brca1_func_score": None,
                    "civic_description": None,
                    "civic_id": None,
                    "civic_score": None,
                    "cosmic_tier": None,
                    "dbsnp_rs_id": None,
                    "dna_change": None,
                    "dp": 20,
                    "ensembl_feature_id": None,
                    "ensembl_gene_id": None,
                    "exon": None,
                    "gc_content": None,
                    "gnomade_af": None,
                    "gnomadg_af": None,
                    "gt": None,
                    "gtex_gene": None,
                    "gwas_disease": None,
                    "gwas_pmid": None,
                    "hess_driver": False,
                    "hess_signature": None,
                    "hgnc_family": None,
                    "hgnc_name": None,
                    "hugo_symbol": None,
                    "intron": None,
                    "lof": None,
                    "molecular_consequence": None,
                    "nmd": None,
                    "oncogene_high_impact": False,
                    "oncokb_effect": "Loss-of-function",
                    "oncokb_hotspot": False,
                    "oncokb_oncogenic": None,
                    "pharmgkb_id": None,
                    "polyphen": None,
                    "protein_change": None,
                    "provean_prediction": None,
                    "ps": None,
                    "ref_count": None,
                    "rescue": True,
                    "revel_score": None,
                    "sift": None,
                    "transcript_likely_lof": None,
                    "tumor_suppressor_high_impact": False,
                    "uniprot_id": None,
                    "variant_info": None,
                    "variant_type": None,
                    "vep_biotype": None,
                    "vep_clin_sig": None,
                    "vep_ensp": None,
                    "vep_existing_variation": None,
                    "vep_hgnc_id": None,
                    "vep_impact": None,
                    "vep_loftool": None,
                    "vep_mane_select": None,
                    "vep_pli_gene_value": None,
                    "vep_somatic": None,
                    "vep_swissprot": None,
                },
            ]
        ).astype(maf_dtypes)

        assert_frame_equal(observed, expected)

    def test_vals(self, db: duckdb.DuckDBPyConnection) -> None:
        db.sql("""
            INSERT INTO variants(
                vid, chrom, pos, ref, alt
            ) VALUES (
                1, 'chr1', 100, 'G', 'C'
            );
            
            INSERT INTO vals_info(
                vid, kind, k, v_varchar, v_integer, v_float, v_integer_arr
            ) VALUES (
                1, 'val', 'af', NULL, NULL, 0.3, NULL
            ), (
                1, 'val', 'dp', NULL, 20, NULL, NULL
            ), (
                1, 'val', 'ad', NULL, 20, NULL, [5, 15]
            ), (
                1, 'val', 'gt', '1|1', NULL, NULL, NULL
            ), (
                1, 'val', 'ps', NULL, 175160788, NULL, NULL
            );
            
            --any rescue state will work for this test
            INSERT INTO vals_info(
                vid, kind, k, v_varchar
            ) VALUES (
                1, 'info', 'oncokb_muteff', 'Loss-of-function'
            );
        """)

        make_views(db)

        observed = get_somatic_variants_as_df(db)

        expected = pd.DataFrame(
            [
                {
                    "chrom": "chr1",
                    "pos": 100,
                    "ref": "G",
                    "alt": "C",
                    "af": 0.3,
                    "alt_count": 15,
                    "am_class": None,
                    "am_pathogenicity": None,
                    "brca1_func_score": None,
                    "civic_description": None,
                    "civic_id": None,
                    "civic_score": None,
                    "cosmic_tier": None,
                    "dbsnp_rs_id": None,
                    "dna_change": None,
                    "dp": 20,
                    "ensembl_feature_id": None,
                    "ensembl_gene_id": None,
                    "exon": None,
                    "gc_content": None,
                    "gnomade_af": None,
                    "gnomadg_af": None,
                    "gt": "1|1",
                    "gtex_gene": None,
                    "gwas_disease": None,
                    "gwas_pmid": None,
                    "hess_driver": False,
                    "hess_signature": None,
                    "hgnc_family": None,
                    "hgnc_name": None,
                    "hugo_symbol": None,
                    "intron": None,
                    "lof": None,
                    "molecular_consequence": None,
                    "nmd": None,
                    "oncogene_high_impact": False,
                    "oncokb_effect": "Loss-of-function",
                    "oncokb_hotspot": False,
                    "oncokb_oncogenic": None,
                    "pharmgkb_id": None,
                    "polyphen": None,
                    "protein_change": None,
                    "provean_prediction": None,
                    "ps": 175160788,
                    "ref_count": 5,
                    "rescue": True,
                    "revel_score": None,
                    "sift": None,
                    "transcript_likely_lof": None,
                    "tumor_suppressor_high_impact": False,
                    "uniprot_id": None,
                    "variant_info": None,
                    "variant_type": None,
                    "vep_biotype": None,
                    "vep_clin_sig": None,
                    "vep_ensp": None,
                    "vep_existing_variation": None,
                    "vep_hgnc_id": None,
                    "vep_impact": None,
                    "vep_loftool": None,
                    "vep_mane_select": None,
                    "vep_pli_gene_value": None,
                    "vep_somatic": None,
                    "vep_swissprot": None,
                },
            ]
        ).astype(maf_dtypes)

        assert_frame_equal(observed, expected)

    def test_info(self, db: duckdb.DuckDBPyConnection) -> None:
        db.sql("""
            INSERT INTO variants(
                vid, chrom, pos, ref, alt
            ) VALUES (
                1, 'chr1', 100, 'G', 'C'
            );
            
            INSERT INTO vals_info(
                vid, kind, k, v_varchar, v_integer, v_float, v_integer_arr
            ) VALUES (
                1, 'val', 'af', NULL, NULL, 0.3, NULL
            ), (
                1, 'val', 'dp', NULL, 20, NULL, NULL
            );
            
            --any rescue state will work for this test
            INSERT INTO vals_info(
                vid, kind, k, v_varchar
            ) VALUES (
                1, 'info', 'oncokb_muteff', 'Loss-of-function'
            );
            
            INSERT INTO vals_info(
                vid, kind, k, v_boolean
            ) VALUES (
                1, 'info', 'oncokb_hotspot', TRUE
            );
            
            INSERT INTO vals_info(
                vid, kind, k, v_varchar
            ) VALUES (
                1, 'info', 'oncokb_oncogenic', 'Oncogenic'
            );

            INSERT INTO vals_info(
                vid, kind, k, v_varchar, v_integer, v_float
            ) VALUES (
                1, 'info', 'civic_desc', 'civic description', NULL, NULL
            ), (
                1, 'info', 'civic_id', NULL, 123, NULL
            ), (
                1, 'info', 'civic_score', NULL, NULL, 45.0
            );

            INSERT INTO vals_info(
                vid, kind, k, v_integer
            ) VALUES (
                1, 'info', 'cmc_tier', 1
            );

            INSERT INTO vals_info(
                vid, kind, k, v_varchar_arr
            ) VALUES (
                1, 'info', 'rs', ['346575']
            );

            INSERT INTO vals_info(
                vid, kind, k, v_float
            ) VALUES (
                1, 'info', 'gc_prop', 0.321
            );

            INSERT INTO vals_info(
                vid, kind, k, v_varchar_arr
            ) VALUES (
                1, 'info', 'mc', ['SO:0001627|intron', 'SO:0001819|synonymous']
            );

            INSERT INTO vals_info(
                vid, kind, k, v_varchar
            ) VALUES (
                1, 'info', 'oc_gtex_gtex_gene', 'LINC00339|LINC00339'
            ), (
                1, 'info', 'oc_gwas_catalog_disease', 'Alzehimer''s'
            ), (
                1, 'info', 'oc_gwas_catalog_pmid', '27863252'
            ), (
                1, 'info', 'oc_pharmgkb_id', 'PA166153763'
            ), (
                1, 'info', 'oc_provean_prediction', 'Damaging'
            );
            
            INSERT INTO vals_info(
                vid, kind, k, v_float
            ) VALUES (
                1, 'info', 'oc_brca1_func_assay_score', -1.5
            ), (
                1, 'info', 'oc_revel_score', 0.012
            );

            INSERT INTO vals_info(
                vid, kind, k, v_varchar
            ) VALUES (
                1, 'info', 'hess', 'sig'
            );
        """)

        make_views(db)

        observed = get_somatic_variants_as_df(db)

        expected = pd.DataFrame(
            [
                {
                    "chrom": "chr1",
                    "pos": 100,
                    "ref": "G",
                    "alt": "C",
                    "af": 0.3,
                    "alt_count": None,
                    "am_class": None,
                    "am_pathogenicity": None,
                    "brca1_func_score": -1.5,
                    "civic_description": "civic description",
                    "civic_id": 123,
                    "civic_score": 45.0,
                    "cosmic_tier": 1,
                    "dbsnp_rs_id": "346575",
                    "dna_change": None,
                    "dp": 20,
                    "ensembl_feature_id": None,
                    "ensembl_gene_id": None,
                    "exon": None,
                    "gc_content": 0.321,
                    "gnomade_af": None,
                    "gnomadg_af": None,
                    "gt": None,
                    "gtex_gene": "LINC00339|LINC00339",
                    "gwas_disease": "Alzehimer's",
                    "gwas_pmid": "27863252",
                    "hess_driver": True,
                    "hess_signature": "sig",
                    "hgnc_family": None,
                    "hgnc_name": None,
                    "hugo_symbol": None,
                    "intron": None,
                    "lof": None,
                    "molecular_consequence": "SO:0001627|intron,SO:0001819|synonymous",
                    "nmd": None,
                    "oncogene_high_impact": False,
                    "oncokb_effect": "Loss-of-function",
                    "oncokb_hotspot": True,
                    "oncokb_oncogenic": "Oncogenic",
                    "pharmgkb_id": "PA166153763",
                    "polyphen": None,
                    "protein_change": None,
                    "provean_prediction": "Damaging",
                    "ps": None,
                    "ref_count": None,
                    "rescue": True,
                    "revel_score": 0.012,
                    "sift": None,
                    "transcript_likely_lof": None,
                    "tumor_suppressor_high_impact": False,
                    "uniprot_id": None,
                    "variant_info": None,
                    "variant_type": None,
                    "vep_biotype": None,
                    "vep_clin_sig": None,
                    "vep_ensp": None,
                    "vep_existing_variation": None,
                    "vep_hgnc_id": None,
                    "vep_impact": None,
                    "vep_loftool": None,
                    "vep_mane_select": None,
                    "vep_pli_gene_value": None,
                    "vep_somatic": None,
                    "vep_swissprot": None,
                },
            ]
        ).astype(maf_dtypes)

        assert_frame_equal(observed, expected)

    def test_vep(self, db: duckdb.DuckDBPyConnection) -> None:
        db.sql("""
            INSERT INTO variants(
                vid, chrom, pos, ref, alt
            ) VALUES (
                1, 'chr1', 100, 'G', 'C'
            );
            
            INSERT INTO vals_info(
                vid, kind, k, v_varchar, v_integer, v_float, v_integer_arr
            ) VALUES (
                1, 'val', 'af', NULL, NULL, 0.3, NULL
            ), (
                1, 'val', 'dp', NULL, 20, NULL, NULL
            );
            
            --any rescue state will work for this test
            INSERT INTO vals_info(
                vid, kind, k, v_varchar
            ) VALUES (
                1, 'info', 'oncokb_muteff', 'Loss-of-function'
            );
            
            INSERT INTO vals_info(
                vid, kind, k, v_boolean, v_json_arr
            ) VALUES (
                1, 'info', 'csq', NULL, ['{
                    "am_class": "likely_benign",
                    "am_pathogenicity": 0.1949,
                    "biotype": "lncRNA",
                    "clin_sig": "benign",
                    "consequence": "intergenic_variant",
                    "ensp": "ENSP00000506999",
                    "existing_variation": "rs10917833",
                    "exon": "1/1",
                    "feature": "ENST00000412705",
                    "gene": "ENSG00000228289",
                    "gnom_ade_af": 0.123,
                    "gnom_adg_af": 0.456,
                    "hgnc_id": "HGNC:3501",
                    "hgvsc": "hgvsc",
                    "hgvsp": "ENSP00000510597.1:p.Asn326Thr",
                    "impact": "HIGH",
                    "intron": "3/4",
                    "loftool": 1.9,
                    "mane_select": "NM_001350197.2",
                    "pli_gene_value": 9.8,
                    "poly_phen": "benign(0)",
                    "sift": "deleterious(0)",
                    "somatic": "0&1",
                    "swissprot": "O60447.171",
                    "symbol": "EVI5",
                    "uniprot_isoform": "O60447-2",
                    "variant_class": "deletion"
                }']
            );
        """)

        make_views(db)

        observed = get_somatic_variants_as_df(db)

        expected = pd.DataFrame(
            [
                {
                    "chrom": "chr1",
                    "pos": 100,
                    "ref": "G",
                    "alt": "C",
                    "af": 0.3,
                    "alt_count": None,
                    "am_class": "likely_benign",
                    "am_pathogenicity": 0.1949,
                    "brca1_func_score": None,
                    "civic_description": None,
                    "civic_id": None,
                    "civic_score": None,
                    "cosmic_tier": None,
                    "dbsnp_rs_id": None,
                    "dna_change": "hgvsc",
                    "dp": 20,
                    "ensembl_feature_id": "ENST00000412705",
                    "ensembl_gene_id": "ENSG00000228289",
                    "exon": "1/1",
                    "gc_content": None,
                    "gnomade_af": 0.123,
                    "gnomadg_af": 0.456,
                    "gt": None,
                    "gtex_gene": None,
                    "gwas_disease": None,
                    "gwas_pmid": None,
                    "hess_driver": False,
                    "hess_signature": None,
                    "hgnc_family": None,
                    "hgnc_name": None,
                    "hugo_symbol": "EVI5",
                    "intron": "3/4",
                    "lof": None,
                    "molecular_consequence": None,
                    "nmd": None,
                    "oncogene_high_impact": False,
                    "oncokb_effect": "Loss-of-function",
                    "oncokb_hotspot": False,
                    "oncokb_oncogenic": None,
                    "pharmgkb_id": None,
                    "polyphen": "benign(0)",
                    "protein_change": "ENSP00000510597.1:p.Asn326Thr",
                    "provean_prediction": None,
                    "ps": None,
                    "ref_count": None,
                    "rescue": True,
                    "revel_score": None,
                    "sift": "deleterious(0)",
                    "transcript_likely_lof": None,
                    "tumor_suppressor_high_impact": False,
                    "uniprot_id": "O60447-2",
                    "variant_info": "intergenic_variant",
                    "variant_type": "deletion",
                    "vep_biotype": "lncRNA",
                    "vep_clin_sig": "benign",
                    "vep_ensp": "ENSP00000506999",
                    "vep_existing_variation": "rs10917833",
                    "vep_hgnc_id": "HGNC:3501",
                    "vep_impact": "HIGH",
                    "vep_loftool": 1.9,
                    "vep_mane_select": "NM_001350197.2",
                    "vep_pli_gene_value": 9.8,
                    "vep_somatic": "0&1",
                    "vep_swissprot": "O60447.171",
                },
            ]
        ).astype(maf_dtypes)

        assert_frame_equal(observed, expected)

    def test_transcript_likely_lof(self, db: duckdb.DuckDBPyConnection) -> None:
        db.sql("""
            INSERT INTO variants(
                vid, chrom, pos, ref, alt
            ) VALUES (
                1, 'chr1', 100, 'G', 'C'
            );
            
            INSERT INTO vals_info(
                vid, kind, k, v_varchar, v_integer, v_float, v_integer_arr
            ) VALUES (
                1, 'val', 'af', NULL, NULL, 0.3, NULL
            ), (
                1, 'val', 'dp', NULL, 20, NULL, NULL
            );
            
            --any rescue state will work for this test
            INSERT INTO vals_info(
                vid, kind, k, v_varchar
            ) VALUES (
                1, 'info', 'oncokb_muteff', 'Loss-of-function'
            );
            
            INSERT INTO vals_info(
                vid, kind, k, v_varchar
            ) VALUES (
                1, 'info', 'oc_revel_all', '[["E1",0.7,0.2],["E2",0.8,0.4],["E3",0.5,0.9]]'
            );
        """)

        make_views(db)

        observed = get_somatic_variants_as_df(db)

        expected = pd.DataFrame(
            [
                {
                    "chrom": "chr1",
                    "pos": 100,
                    "ref": "G",
                    "alt": "C",
                    "af": 0.3,
                    "alt_count": None,
                    "am_class": None,
                    "am_pathogenicity": None,
                    "brca1_func_score": None,
                    "civic_description": None,
                    "civic_id": None,
                    "civic_score": None,
                    "cosmic_tier": None,
                    "dbsnp_rs_id": None,
                    "dna_change": None,
                    "dp": 20,
                    "ensembl_feature_id": None,
                    "ensembl_gene_id": None,
                    "exon": None,
                    "gc_content": None,
                    "gnomade_af": None,
                    "gnomadg_af": None,
                    "gt": None,
                    "gtex_gene": None,
                    "gwas_disease": None,
                    "gwas_pmid": None,
                    "hess_driver": False,
                    "hess_signature": None,
                    "hgnc_family": None,
                    "hgnc_name": None,
                    "hugo_symbol": None,
                    "intron": None,
                    "lof": None,
                    "molecular_consequence": None,
                    "nmd": None,
                    "oncogene_high_impact": False,
                    "oncokb_effect": "Loss-of-function",
                    "oncokb_hotspot": False,
                    "oncokb_oncogenic": None,
                    "pharmgkb_id": None,
                    "polyphen": None,
                    "protein_change": None,
                    "provean_prediction": None,
                    "ps": None,
                    "ref_count": None,
                    "rescue": True,
                    "revel_score": None,
                    "sift": None,
                    "transcript_likely_lof": "E1;E2",
                    "tumor_suppressor_high_impact": False,
                    "uniprot_id": None,
                    "variant_info": None,
                    "variant_type": None,
                    "vep_biotype": None,
                    "vep_clin_sig": None,
                    "vep_ensp": None,
                    "vep_existing_variation": None,
                    "vep_hgnc_id": None,
                    "vep_impact": None,
                    "vep_loftool": None,
                    "vep_mane_select": None,
                    "vep_pli_gene_value": None,
                    "vep_somatic": None,
                    "vep_swissprot": None,
                },
            ]
        ).astype(maf_dtypes)

        assert_frame_equal(observed, expected)

    def test_nmd(self, db: duckdb.DuckDBPyConnection) -> None:
        db.sql("""
            INSERT INTO variants(
                vid, chrom, pos, ref, alt
            ) VALUES (
                1, 'chr1', 100, 'G', 'C'
            );
            
            INSERT INTO vals_info(
                vid, kind, k, v_varchar, v_integer, v_float, v_integer_arr
            ) VALUES (
                1, 'val', 'af', NULL, NULL, 0.3, NULL
            ), (
                1, 'val', 'dp', NULL, 20, NULL, NULL
            );
            
            --any rescue state will work for this test
            INSERT INTO vals_info(
                vid, kind, k, v_varchar
            ) VALUES (
                1, 'info', 'oncokb_muteff', 'Loss-of-function'
            );
            
            INSERT INTO vals_info(
                vid, kind, k, v_json_arr
            ) VALUES (
                1, 'info', 'nmd', ['{"gene_name":"LRRC20","gene_id":"E1.15","number_of_transcripts_in_gene":"1","percent_of_transcripts_affected":"1.00"}']
            );
        """)

        make_views(db)

        observed = get_somatic_variants_as_df(db)

        expected = pd.DataFrame(
            [
                {
                    "chrom": "chr1",
                    "pos": 100,
                    "ref": "G",
                    "alt": "C",
                    "af": 0.3,
                    "alt_count": None,
                    "am_class": None,
                    "am_pathogenicity": None,
                    "brca1_func_score": None,
                    "civic_description": None,
                    "civic_id": None,
                    "civic_score": None,
                    "cosmic_tier": None,
                    "dbsnp_rs_id": None,
                    "dna_change": None,
                    "dp": 20,
                    "ensembl_feature_id": None,
                    "ensembl_gene_id": None,
                    "exon": None,
                    "gc_content": None,
                    "gnomade_af": None,
                    "gnomadg_af": None,
                    "gt": None,
                    "gtex_gene": None,
                    "gwas_disease": None,
                    "gwas_pmid": None,
                    "hess_driver": False,
                    "hess_signature": None,
                    "hgnc_family": None,
                    "hgnc_name": None,
                    "hugo_symbol": None,
                    "intron": None,
                    "lof": None,
                    "molecular_consequence": None,
                    "nmd": """[{"gene_name":"LRRC20","gene_id":"E1.15","number_of_transcripts_in_gene":"1","percent_of_transcripts_affected":"1.00"}]""",
                    "oncogene_high_impact": False,
                    "oncokb_effect": "Loss-of-function",
                    "oncokb_hotspot": False,
                    "oncokb_oncogenic": None,
                    "pharmgkb_id": None,
                    "polyphen": None,
                    "protein_change": None,
                    "provean_prediction": None,
                    "ps": None,
                    "ref_count": None,
                    "rescue": True,
                    "revel_score": None,
                    "sift": None,
                    "transcript_likely_lof": None,
                    "tumor_suppressor_high_impact": False,
                    "uniprot_id": None,
                    "variant_info": None,
                    "variant_type": None,
                    "vep_biotype": None,
                    "vep_clin_sig": None,
                    "vep_ensp": None,
                    "vep_existing_variation": None,
                    "vep_hgnc_id": None,
                    "vep_impact": None,
                    "vep_loftool": None,
                    "vep_mane_select": None,
                    "vep_pli_gene_value": None,
                    "vep_somatic": None,
                    "vep_swissprot": None,
                },
            ]
        ).astype(maf_dtypes)

        assert_frame_equal(observed, expected)

    def test_hgnc(self, db: duckdb.DuckDBPyConnection) -> None:
        db.sql("""
            INSERT INTO variants(
                vid, chrom, pos, ref, alt
            ) VALUES (
                1, 'chr1', 100, 'G', 'C'
            ), (
                2, 'chr2', 200, 'A', 'T'
            );
            
            INSERT INTO vals_info(
                vid, kind, k, v_varchar, v_integer, v_float, v_integer_arr
            ) VALUES (
                1, 'val', 'af', NULL, NULL, 0.3, NULL
            ), (
                1, 'val', 'dp', NULL, 20, NULL, NULL
            ), (
                2, 'val', 'af', NULL, NULL, 0.3, NULL
            ), (
                2, 'val', 'dp', NULL, 20, NULL, NULL
            );
            
            --any rescue state will work for this test
            INSERT INTO vals_info(
                vid, kind, k, v_varchar
            ) VALUES (
                1, 'info', 'oncokb_muteff', 'Loss-of-function'
            ), (
                2, 'info', 'oncokb_muteff', 'Loss-of-function'
            );
            
            INSERT INTO vals_info(
                vid, kind, k, v_varchar_arr
            ) VALUES (
                1, 'info', 'hgnc_name', ['name1', 'name2']
            ), (
                1, 'info', 'hgnc_group', ['group1']
            ), (
                2, 'info', 'hgnc_name', ['nameA']
            );
        """)

        make_views(db)

        observed = get_somatic_variants_as_df(db)

        expected = pd.DataFrame(
            [
                {
                    "chrom": "chr1",
                    "pos": 100,
                    "ref": "G",
                    "alt": "C",
                    "af": 0.3,
                    "alt_count": None,
                    "am_class": None,
                    "am_pathogenicity": None,
                    "brca1_func_score": None,
                    "civic_description": None,
                    "civic_id": None,
                    "civic_score": None,
                    "cosmic_tier": None,
                    "dbsnp_rs_id": None,
                    "dna_change": None,
                    "dp": 20,
                    "ensembl_feature_id": None,
                    "ensembl_gene_id": None,
                    "exon": None,
                    "gc_content": None,
                    "gnomade_af": None,
                    "gnomadg_af": None,
                    "gt": None,
                    "gtex_gene": None,
                    "gwas_disease": None,
                    "gwas_pmid": None,
                    "hess_driver": False,
                    "hess_signature": None,
                    "hgnc_family": "group1",
                    "hgnc_name": "name1;name2",
                    "hugo_symbol": None,
                    "intron": None,
                    "lof": None,
                    "molecular_consequence": None,
                    "nmd": None,
                    "oncogene_high_impact": False,
                    "oncokb_effect": "Loss-of-function",
                    "oncokb_hotspot": False,
                    "oncokb_oncogenic": None,
                    "pharmgkb_id": None,
                    "polyphen": None,
                    "protein_change": None,
                    "provean_prediction": None,
                    "ps": None,
                    "ref_count": None,
                    "rescue": True,
                    "revel_score": None,
                    "sift": None,
                    "transcript_likely_lof": None,
                    "tumor_suppressor_high_impact": False,
                    "uniprot_id": None,
                    "variant_info": None,
                    "variant_type": None,
                    "vep_biotype": None,
                    "vep_clin_sig": None,
                    "vep_ensp": None,
                    "vep_existing_variation": None,
                    "vep_hgnc_id": None,
                    "vep_impact": None,
                    "vep_loftool": None,
                    "vep_mane_select": None,
                    "vep_pli_gene_value": None,
                    "vep_somatic": None,
                    "vep_swissprot": None,
                },
                {
                    "chrom": "chr2",
                    "pos": 200,
                    "ref": "A",
                    "alt": "T",
                    "af": 0.3,
                    "alt_count": None,
                    "am_class": None,
                    "am_pathogenicity": None,
                    "brca1_func_score": None,
                    "civic_description": None,
                    "civic_id": None,
                    "civic_score": None,
                    "cosmic_tier": None,
                    "dbsnp_rs_id": None,
                    "dna_change": None,
                    "dp": 20,
                    "ensembl_feature_id": None,
                    "ensembl_gene_id": None,
                    "exon": None,
                    "gc_content": None,
                    "gnomade_af": None,
                    "gnomadg_af": None,
                    "gt": None,
                    "gtex_gene": None,
                    "gwas_disease": None,
                    "gwas_pmid": None,
                    "hess_driver": False,
                    "hess_signature": None,
                    "hgnc_family": None,
                    "hgnc_name": "nameA",
                    "hugo_symbol": None,
                    "intron": None,
                    "lof": None,
                    "molecular_consequence": None,
                    "nmd": None,
                    "oncogene_high_impact": False,
                    "oncokb_effect": "Loss-of-function",
                    "oncokb_hotspot": False,
                    "oncokb_oncogenic": None,
                    "pharmgkb_id": None,
                    "polyphen": None,
                    "protein_change": None,
                    "provean_prediction": None,
                    "ps": None,
                    "ref_count": None,
                    "rescue": True,
                    "revel_score": None,
                    "sift": None,
                    "transcript_likely_lof": None,
                    "tumor_suppressor_high_impact": False,
                    "uniprot_id": None,
                    "variant_info": None,
                    "variant_type": None,
                    "vep_biotype": None,
                    "vep_clin_sig": None,
                    "vep_ensp": None,
                    "vep_existing_variation": None,
                    "vep_hgnc_id": None,
                    "vep_impact": None,
                    "vep_loftool": None,
                    "vep_mane_select": None,
                    "vep_pli_gene_value": None,
                    "vep_somatic": None,
                    "vep_swissprot": None,
                },
            ]
        ).astype(maf_dtypes)

        assert_frame_equal(observed, expected)
