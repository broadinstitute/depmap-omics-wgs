from tempfile import NamedTemporaryFile, TemporaryDirectory

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
def db_setup() -> None:
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
                
                CREATE TABLE vals(
                    vid UINTEGER,
                    k VARCHAR NOT NULL,
                    v_varchar VARCHAR,
                    v_integer INTEGER,
                    v_float FLOAT,
                    v_integer_arr INTEGER[],
                    FOREIGN KEY (vid) REFERENCES variants(vid)
                );
                
                CREATE TABLE info(
                    vid UINTEGER,
                    k VARCHAR NOT NULL,
                    v_varchar VARCHAR,
                    v_integer INTEGER,
                    v_float FLOAT,
                    v_boolean BOOLEAN,
                    v_varchar_arr VARCHAR[],
                    v_integer_arr INTEGER[],
                    v_json_arr JSON[],
                    FOREIGN KEY (vid) REFERENCES variants(vid)
                );
            """)

            yield db


@pytest.fixture
def db(db_setup: duckdb.DuckDBPyConnection) -> None:
    db_setup.sql("""
        DROP TABLE IF EXISTS somatic_variants;
        DROP TABLE IF EXISTS vep;
        
        TRUNCATE TABLE info;
        TRUNCATE TABLE vals;
        TRUNCATE TABLE variants;
    """)

    yield db_setup


@pytest.mark.usefixtures("db_setup")
class Test:
    def test_noop(self, db: duckdb.DuckDBPyConnection) -> None:
        make_views(db)

        observed = get_somatic_variants_as_df(db)
        expected = pd.DataFrame([], columns=list(maf_dtypes.keys())).astype(maf_dtypes)
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
            
            INSERT INTO vals(
                vid, k, v_integer, v_float
            ) VALUES (
                1, 'af', NULL, 0.3
            ), (
                1, 'dp', 20, NULL
            ), (
                2, 'af', NULL, 0.3
            ), (
                2, 'dp', 20, NULL
            ), (
                3, 'af', NULL, 0.3
            ), (
                3, 'dp', 20, NULL
            );
            
            --any rescue state will work for this test
            INSERT INTO info(
                vid, k, v_varchar
            ) VALUES (
                1, 'oncokb_muteff', 'Loss-of-function'
            ), (
                1, 'as_filter_status', 'SITE|SITE'
            ), (
                2, 'oncokb_muteff', 'Loss-of-function'
            ), (
                3, 'oncokb_muteff', 'Loss-of-function'
            ), (
                3, 'as_filter_status', 'strand_bias|SITE'
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
                    "lof_gene_id": None,
                    "lof_gene_name": None,
                    "lof_number_of_transcripts_in_gene": None,
                    "lof_prop_of_transcripts_affected": None,
                    "molecular_consequence": None,
                    "nmd": None,
                    "oncogene_high_impact": False,
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
            
            INSERT INTO vals(
                vid, k, v_integer, v_float
            ) VALUES (
                1, 'af', NULL, 0.3
            ), (
                1, 'dp', 20, NULL
            ), (
                2, 'af', NULL, 0.3
            ), (
                2, 'dp', 20, NULL
            ), (
                3, 'af', NULL, 0.3
            ), (
                3, 'dp', 20, NULL
            );
            
            INSERT INTO info(
                vid, k, v_varchar
            ) VALUES (
                1, 'oncokb_muteff', 'Loss-of-function'
            ), (
                2, 'oncokb_muteff', 'Gain-of-function'
            ), (
                3, 'oncokb_muteff', 'other'
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
                    "lof_gene_id": None,
                    "lof_gene_name": None,
                    "lof_number_of_transcripts_in_gene": None,
                    "lof_prop_of_transcripts_affected": None,
                    "molecular_consequence": None,
                    "nmd": None,
                    "oncogene_high_impact": False,
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
                    "lof_gene_id": None,
                    "lof_gene_name": None,
                    "lof_number_of_transcripts_in_gene": None,
                    "lof_prop_of_transcripts_affected": None,
                    "molecular_consequence": None,
                    "nmd": None,
                    "oncogene_high_impact": False,
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
            
            INSERT INTO vals(
                vid, k, v_integer, v_float
            ) VALUES (
                1, 'af', NULL, 0.3
            ), (
                1, 'dp', 20, NULL
            ), (
                2, 'af', NULL, 0.3
            ), (
                2, 'dp', 20, NULL
            );
            
            INSERT INTO info(
                vid, k, v_varchar
            ) VALUES (
                1, 'oncokb_oncogenic', 'Oncogenic'
            ), (
                2, 'oncokb_oncogenic', 'not'
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
                    "lof_gene_id": None,
                    "lof_gene_name": None,
                    "lof_number_of_transcripts_in_gene": None,
                    "lof_prop_of_transcripts_affected": None,
                    "molecular_consequence": None,
                    "nmd": None,
                    "oncogene_high_impact": False,
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
            
            INSERT INTO vals(
                vid, k, v_integer, v_float
            ) VALUES (
                1, 'af', NULL, 0.3
            ), (
                1, 'dp', 20, NULL
            ), (
                2, 'af', NULL, 0.3
            ), (
                2, 'dp', 20, NULL
            );
            
            INSERT INTO info(
                vid, k, v_boolean
            ) VALUES (
                1, 'oncokb_hotspot', TRUE
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
                    "lof_gene_id": None,
                    "lof_gene_name": None,
                    "lof_number_of_transcripts_in_gene": None,
                    "lof_prop_of_transcripts_affected": None,
                    "molecular_consequence": None,
                    "nmd": None,
                    "oncogene_high_impact": False,
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
            
            INSERT INTO vals(
                vid, k, v_integer, v_float
            ) VALUES (
                1, 'af', NULL, 0.3
            ), (
                1, 'dp', 20, NULL
            ), (
                2, 'af', NULL, 0.3
            ), (
                2, 'dp', 20, NULL
            );
            
            INSERT INTO info(
                vid, k, v_integer
            ) VALUES (
                1, 'cmc_tier', 1 
            ), (
                2, 'cmc_tier', 3 
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
                    "lof_gene_id": None,
                    "lof_gene_name": None,
                    "lof_number_of_transcripts_in_gene": None,
                    "lof_prop_of_transcripts_affected": None,
                    "molecular_consequence": None,
                    "nmd": None,
                    "oncogene_high_impact": False,
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
            
            INSERT INTO vals(
                vid, k, v_integer, v_float
            ) VALUES (
                1, 'af', NULL, 0.3
            ), (
                1, 'dp', 20, NULL
            ), (
                2, 'af', NULL, 0.3
            ), (
                2, 'dp', 20, NULL
            );
            
            INSERT INTO info(
                vid, k, v_float
            ) VALUES (
                1, 'oc_brca1_func_assay_score', -1.5
            ), (
                2, 'oc_brca1_func_assay_score', -0.9
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
                    "lof_gene_id": None,
                    "lof_gene_name": None,
                    "lof_number_of_transcripts_in_gene": None,
                    "lof_prop_of_transcripts_affected": None,
                    "molecular_consequence": None,
                    "nmd": None,
                    "oncogene_high_impact": False,
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
            
            INSERT INTO vals(
                vid, k, v_integer, v_float
            ) VALUES (
                1, 'af', NULL, 0.3
            ), (
                1, 'dp', 20, NULL
            ), (
                2, 'af', NULL, 0.3
            ), (
                2, 'dp', 20, NULL
            ), (
                3, 'af', NULL, 0.3
            ), (
                3, 'dp', 20, NULL
            );
            
            INSERT INTO info(
                vid, k, v_boolean, v_json_arr
            ) VALUES (
                1, 'oncogene', TRUE, NULL 
            ), (
                1, 'csq', NULL, ['{"impact": "HIGH"}']
            ), (
                2, 'oncogene', TRUE, NULL 
            ), (
                2, 'csq', NULL, ['{"impact": "LOW"}']
            ), (
                3, 'oncogene', TRUE, NULL 
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
                    "lof_gene_id": None,
                    "lof_gene_name": None,
                    "lof_number_of_transcripts_in_gene": None,
                    "lof_prop_of_transcripts_affected": None,
                    "molecular_consequence": None,
                    "nmd": None,
                    "oncogene_high_impact": True,
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
            
            INSERT INTO vals(
                vid, k, v_integer, v_float
            ) VALUES (
                1, 'af', NULL, 0.3
            ), (
                1, 'dp', 20, NULL
            ), (
                2, 'af', NULL, 0.3
            ), (
                2, 'dp', 20, NULL
            ), (
                3, 'af', NULL, 0.3
            ), (
                3, 'dp', 20, NULL
            );
            
            INSERT INTO info(
                vid, k, v_boolean, v_json_arr
            ) VALUES (
                1, 'tsg', TRUE, NULL 
            ), (
                1, 'csq', NULL, ['{"impact": "HIGH"}']
            ), (
                2, 'tsg', TRUE, NULL 
            ), (
                2, 'csq', NULL, ['{"impact": "LOW"}']
            ), (
                3, 'tsg', TRUE, NULL 
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
                    "lof_gene_id": None,
                    "lof_gene_name": None,
                    "lof_number_of_transcripts_in_gene": None,
                    "lof_prop_of_transcripts_affected": None,
                    "molecular_consequence": None,
                    "nmd": None,
                    "oncogene_high_impact": False,
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
            
            INSERT INTO vals(
                vid, k, v_integer, v_float
            ) VALUES (
                1, 'af', NULL, 0.3
            ), (
                1, 'dp', 20, NULL
            ), (
                2, 'af', NULL, 0.3
            ), (
                2, 'dp', 20, NULL
            );
            
            INSERT INTO info(
                vid, k, v_varchar
            ) VALUES (
                1, 'hess', 'sig'
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
                    "lof_gene_id": None,
                    "lof_gene_name": None,
                    "lof_number_of_transcripts_in_gene": None,
                    "lof_prop_of_transcripts_affected": None,
                    "molecular_consequence": None,
                    "nmd": None,
                    "oncogene_high_impact": False,
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
            
            INSERT INTO vals(
                vid, k, v_integer, v_float
            ) VALUES (
                1, 'af', NULL, 0.3
            ), (
                1, 'dp', 20, NULL
            ), (
                2, 'af', NULL, 0.3
            ), (
                2, 'dp', 20, NULL
            );
            
            INSERT INTO info(
                vid, k, v_json_arr
            ) VALUES (
                1, 'csq', ['{"symbol": "TERT"}']
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
                    "lof_gene_id": None,
                    "lof_gene_name": None,
                    "lof_number_of_transcripts_in_gene": None,
                    "lof_prop_of_transcripts_affected": None,
                    "molecular_consequence": None,
                    "nmd": None,
                    "oncogene_high_impact": False,
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
            
            INSERT INTO vals(
                vid, k, v_integer, v_float
            ) VALUES (
                1, 'af', NULL, 0.3
            ), (
                1, 'dp', 20, NULL
            ), (
                2, 'af', NULL, 0.3
            ), (
                2, 'dp', 20, NULL
            );
            
            INSERT INTO info(
                vid, k, v_json_arr
            ) VALUES (
                1, 'csq', ['{"symbol": "MET"}']
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
                    "lof_gene_id": None,
                    "lof_gene_name": None,
                    "lof_number_of_transcripts_in_gene": None,
                    "lof_prop_of_transcripts_affected": None,
                    "molecular_consequence": None,
                    "nmd": None,
                    "oncogene_high_impact": False,
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
            
            INSERT INTO vals(
                vid, k, v_integer, v_float
            ) VALUES (
                1, 'af', NULL, 0.3
            ), (
                1, 'dp', 20, NULL
            ), (
                2, 'af', NULL, 0.3
            ), (
                2, 'dp', 20, NULL
            ), (
                3, 'af', NULL, 0.3
            ), (
                3, 'dp', 20, NULL
            );
            
            INSERT INTO info(
                vid, k, v_json_arr
            ) VALUES (
                1, 'csq', ['{"impact": "HIGH", "consequence": "foo&splice_region"}']
            ), (
                2, 'csq', ['{"impact": "HIGH", "consequence": "other"}']
            ), (
                3, 'csq', ['{"impact": "LOW", "consequence": "splice"}']
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
                    "lof_gene_id": None,
                    "lof_gene_name": None,
                    "lof_number_of_transcripts_in_gene": None,
                    "lof_prop_of_transcripts_affected": None,
                    "molecular_consequence": None,
                    "nmd": None,
                    "oncogene_high_impact": False,
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
            
            INSERT INTO vals(
                vid, k, v_integer, v_float
            ) VALUES (
                1, 'af', NULL, 0.3
            ), (
                1, 'dp', 20, NULL
            ), (
                2, 'af', NULL, 0.3
            ), (
                2, 'dp', 20, NULL
            ), (
                3, 'af', NULL, 0.3
            ), (
                3, 'dp', 20, NULL
            );
            
            INSERT INTO info(
                vid, k, v_json_arr
            ) VALUES (
                1, 'csq', ['{"impact": "HIGH", "hgvsp": "ENSP15.3:p.Glu32del"}']
            ), (
                2, 'csq', ['{"impact": "HIGH", "hgvsp": "ENSP30.3:p.Asp637="}']
            ), (
                3, 'csq', ['{"impact": "LOW", "hgvsp": null}']
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
                    "lof_gene_id": None,
                    "lof_gene_name": None,
                    "lof_number_of_transcripts_in_gene": None,
                    "lof_prop_of_transcripts_affected": None,
                    "molecular_consequence": None,
                    "nmd": None,
                    "oncogene_high_impact": False,
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
            
            INSERT INTO vals(
                vid, k, v_integer, v_float
            ) VALUES (
                1, 'af', NULL, 0.3
            ), (
                1, 'dp', 20, NULL
            ), (
                2, 'af', NULL, 0.3
            ), (
                2, 'dp', 20, NULL
            ), (
                3, 'af', NULL, 0.3
            ), (
                3, 'dp', 20, NULL
            );
            
            INSERT INTO info(
                vid, k, v_json_arr
            ) VALUES (
                1, 'csq', ['{"impact": "HIGH", "consequence": "splice"}']
            ), (
                2, 'csq', ['{"impact": "HIGH", "consequence": "splice"}']
            ), (
                3, 'csq', ['{"impact": "HIGH", "consequence": "splice"}']
            );
            
            -- undo some of the previous positive selections
            INSERT INTO info(
                vid, k, v_boolean
            ) VALUES (
                2, 'segdup', TRUE
            ), (
                3, 'rm', TRUE
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
                    "lof_gene_id": None,
                    "lof_gene_name": None,
                    "lof_number_of_transcripts_in_gene": None,
                    "lof_prop_of_transcripts_affected": None,
                    "molecular_consequence": None,
                    "nmd": None,
                    "oncogene_high_impact": False,
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
            
            INSERT INTO vals(
                vid, k, v_integer, v_float
            ) VALUES (
                1, 'af', NULL, 0.3
            ), (
                1, 'dp', 20, NULL
            ), (
                2, 'af', NULL, 0.3
            ), (
                2, 'dp', 20, NULL
            ), (
                3, 'af', NULL, 0.3
            ), (
                3, 'dp', 20, NULL
            ), (
                4, 'af', NULL, 0.3
            ), (
                4, 'dp', 20, NULL
            );
            
            INSERT INTO info(
                vid, k, v_json_arr
            ) VALUES (
                1, 'csq', ['{"impact": "HIGH", "consequence": "splice", "gnom_ade_af": 0.000001, "gnom_adg_af": 0.000002}']
            ), (
                2, 'csq', ['{"impact": "HIGH", "consequence": "splice", "gnom_ade_af": 0.000001, "gnom_adg_af": 0.1}']
            ), (
                3, 'csq', ['{"impact": "HIGH", "consequence": "splice", "gnom_ade_af": null, "gnom_adg_af": 0.1}']
            ), (
                4, 'csq', ['{"impact": "HIGH", "consequence": "splice", "gnom_ade_af": 0.000001, "gnom_adg_af": 0.000002}']
            );
            
            INSERT INTO info(
                vid, k, v_boolean
            ) VALUES (
                4, 'pon', TRUE
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
                    "lof_gene_id": None,
                    "lof_gene_name": None,
                    "lof_number_of_transcripts_in_gene": None,
                    "lof_prop_of_transcripts_affected": None,
                    "molecular_consequence": None,
                    "nmd": None,
                    "oncogene_high_impact": False,
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
