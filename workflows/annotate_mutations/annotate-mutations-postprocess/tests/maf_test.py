import os
from pathlib import Path
from tempfile import NamedTemporaryFile, TemporaryDirectory

import duckdb
import pandas as pd
from pandas.testing import assert_frame_equal

from annotate_mutations_postprocess.maf import convert_duckdb_to_maf

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


class TestMaf:
    def test_noop(self):
        with TemporaryDirectory() as tmpdir:
            db_path = NamedTemporaryFile(dir=tmpdir, suffix="duckdb").name
            parquet_dir_path = os.path.join(tmpdir, "parquet")
            out_file_path = NamedTemporaryFile(dir=tmpdir, suffix="parquet").name

            with duckdb.connect(db_path) as db:
                db.sql("""
                    CREATE TABLE val_info_types(
                        id VARCHAR,
                        number VARCHAR,
                        "type" VARCHAR,
                        description VARCHAR,
                        kind VARCHAR,
                        has_children BOOLEAN,
                        parent_id VARCHAR,
                        ix UBIGINT,
                        col_def VARCHAR,
                        is_list BOOLEAN,
                        v_col_name VARCHAR,
                        id_snake VARCHAR
                    );
                    
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
                    
                    CREATE TABLE vals(vid UINTEGER,
                        k VARCHAR NOT NULL,
                        v_integer_arr INTEGER[],
                        v_float FLOAT,
                        v_integer INTEGER,
                        v_varchar VARCHAR,
                        FOREIGN KEY (vid) REFERENCES variants(vid)
                    );
                    
                    CREATE TABLE info(
                        vid UINTEGER,
                        k VARCHAR NOT NULL,
                        v_varchar VARCHAR,
                        v_integer INTEGER,
                        v_float FLOAT,
                        v_integer_arr INTEGER[],
                        v_boolean BOOLEAN,
                        v_varchar_arr VARCHAR[],
                        v_json_arr JSON[],
                        FOREIGN KEY (vid) REFERENCES variants(vid)
                    );
                """)

                db.sql(f"EXPORT DATABASE '{parquet_dir_path}' (FORMAT PARQUET)")

            os.remove(db_path)

            convert_duckdb_to_maf(
                Path(db_path), Path(parquet_dir_path), Path(out_file_path)
            )

            observed = pd.read_parquet(out_file_path)
            expected = pd.DataFrame([], columns=list(maf_dtypes.keys())).astype(
                maf_dtypes
            )
            assert_frame_equal(observed, expected)
