import logging
import os
import re
import subprocess
from math import ceil
from pathlib import Path
from urllib.parse import unquote

import duckdb
import pandas as pd
from caseconverter import snakecase
from duckdb import DuckDBPyConnection
from vcf_info_merger.merge import get_header_lines


def create_and_populate_db(
    vcf_path: Path,
    db_path: Path,
    compound_info_fields: set[str],
    info_cols_ignored: set[str],
    url_encoded_col_name_regex: re.Pattern,
) -> None:
    val_info_types = get_vcf_val_info_types(
        vcf_path, compound_info_fields, info_cols_ignored
    )

    try:
        os.remove(db_path)
    except OSError:
        pass

    with duckdb.connect(db_path) as db:
        set_up_db(db, val_info_types)

        tab_path = vcf_path.with_suffix("").with_suffix(".tsv")
        write_tab_vcf(vcf_gz_path=vcf_path, tab_path=tab_path)

        logging.info(f"Reading {vcf_path} into {db_path}")
        populate_db(db, tab_path, val_info_types)
        pass


def set_up_db(db: DuckDBPyConnection, val_info_types: pd.DataFrame) -> None:
    db.sql("""
        CREATE TABLE IF NOT EXISTS vcf_lines (
            chrom VARCHAR NOT NULL,
            pos UINTEGER NOT NULL,
            id VARCHAR,
            ref VARCHAR,
            alt VARCHAR,
            qual VARCHAR,
            filters VARCHAR,
            info VARCHAR,
            format VARCHAR,
            values VARCHAR
        );
    """)

    db.sql("""
        CREATE TABLE IF NOT EXISTS variants (
            vid VARCHAR PRIMARY KEY,
            chrom VARCHAR NOT NULL,
            pos UINTEGER NOT NULL,
            id VARCHAR,
            ref VARCHAR,
            alt VARCHAR,
            qual VARCHAR,
            filters VARCHAR[]
        );
    """)

    db.sql("""
        CREATE TABLE IF NOT EXISTS kv (
            vid VARCHAR,
            k VARCHAR,
            v VARCHAR
        );
    """)

    db.sql("""
        CREATE TABLE IF NOT EXISTS kv_compound_info (
            vid VARCHAR,
            k VARCHAR,
            k_sub VARCHAR,
            ix INTEGER,
            v VARCHAR
        );
    """)

    val_types = val_info_types.loc[val_info_types["kind"].eq("value")].copy()

    req_val_col_types = val_types[["v_col_name", "col_def"]].drop_duplicates()
    val_cols = req_val_col_types["v_col_name"] + " " + req_val_col_types["col_def"]

    db.sql(f"""
        CREATE TABLE IF NOT EXISTS vals (
            vid VARCHAR,
            k VARCHAR,
            {', '.join(val_cols)}
        );
    """)

    db.sql(f"""
        CREATE TABLE IF NOT EXISTS vals_fk (
            vid VARCHAR REFERENCES variants (vid),
            k VARCHAR NOT NULL,
            {', '.join(val_cols)}
        );
    """)

    info_types = val_info_types.loc[val_info_types["kind"].eq("info")].copy()

    req_info_col_types = info_types[["v_col_name", "col_def"]].drop_duplicates()
    info_cols = req_info_col_types["v_col_name"] + " " + req_info_col_types["col_def"]

    db.sql(f"""
        CREATE TABLE IF NOT EXISTS info (
            vid VARCHAR,
            k VARCHAR,
            {', '.join(info_cols)}
        );
    """)

    db.sql(f"""
        CREATE TABLE IF NOT EXISTS info_fk (
            vid VARCHAR REFERENCES variants (vid),
            k VARCHAR NOT NULL,
            {', '.join(info_cols)}
        );
    """)

    sub_fields = val_info_types.loc[
        val_info_types["parent_id_db"].isin(info_types["id_db"]),
        ["ix", "parent_id_db", "id_db"],
    ].rename(columns={"parent_id_db": "k"})

    db.register("sub_fields", sub_fields)


def get_vcf_val_info_types(
    path: Path,
    compound_info_fields: set[str],
    info_cols_ignored: set[str],
) -> pd.DataFrame:
    header_lines = get_header_lines([path])
    header_lines = [
        x for x in header_lines if x.startswith("##FORMAT") or x.startswith("##INFO")
    ]

    type_map = {
        "Integer": "INTEGER",
        "Float": "FLOAT",
        "String": "VARCHAR",
        "Character": "VARCHAR",
        "Flag": "BOOLEAN",
    }

    arr = []
    sub_arr = []

    for x in header_lines:
        kind = re.search(r"^##(\w+)", x).group(1).lower()
        interior = re.search(r"<(.+)>$", x).group(1)
        parts = re.findall(r'([A-Za-z0-9_]+)=(".*?"|[^,]+)', interior)

        d = {k.lower(): v.strip('"') for k, v in parts}
        d["kind"] = kind if kind == "info" else "value"

        if d["kind"] == "info" and d["id"] in info_cols_ignored:
            continue

        d["id_db"] = d["id"]

        if d["id"] in compound_info_fields:
            d["has_children"] = True
            d["type"] = "JSON"
            # d["type"] = f"info_{d['id_db']}"

            desc = re.search(r"^.+:['\s]*([^']+)['\s]*$", d["description"]).group(1)
            subfields = re.split(r"\s*\|\s*", desc)

            for ix, s in enumerate(subfields):
                dsub = {
                    "id": s,
                    "id_db": s,
                    "has_children": False,
                    "number": "1",
                    "type": "VARCHAR",
                    "kind": "sub_info",
                    "parent_id": d["id"],
                    "parent_id_db": d["id_db"],
                    "ix": ix + 1,
                }

                sub_arr.append(dsub)

        else:
            d["has_children"] = False
            d["type"] = type_map[d["type"]]

        arr.append(d)

    df = pd.DataFrame(arr + sub_arr).astype(
        {
            "id": "string",
            "number": "string",
            "type": "string",
            "description": "string",
            "kind": "string",
            "id_db": "string",
            "has_children": "boolean",
            "parent_id": "string",
            "parent_id_db": "string",
            "ix": "Int64",
        }
    )

    df["col_def"] = df["type"]
    df["is_list"] = ~df["number"].isin({"0", "1"})
    df.loc[df["is_list"], "col_def"] += "[]"

    df["v_col_name"] = "v_" + df["type"].str.lower()
    df.loc[df["is_list"], "v_col_name"] += "_arr"

    return df


def write_tab_vcf(vcf_gz_path: Path, tab_path: Path) -> None:
    logging.info(f"Converting {vcf_gz_path} to TSV")
    subprocess.run(
        ["bcftools", "view", vcf_gz_path, "--no-header", "-o", tab_path]
        # ["bcftools", "view", vcf_gz_path, "-r", "chr1", "--no-header", "-o", tab_path]
    )


def populate_db(
    db: DuckDBPyConnection, tab_path: Path, val_info_types: pd.DataFrame
) -> None:
    db.sql(f"""
        COPY
            vcf_lines
        FROM
            '{tab_path}' (DELIMITER '\\t', AUTO_DETECT false);
    """)

    db.sql("""
        UPDATE vcf_lines SET id = NULL where id = '.';
        UPDATE vcf_lines SET ref = NULL where ref = '.';
        UPDATE vcf_lines SET alt = NULL where alt = '.';
        UPDATE vcf_lines SET qual = NULL where qual = '.';
        UPDATE vcf_lines SET filters = NULL where filters = '.';
        UPDATE vcf_lines SET info = NULL where info = '.';
        UPDATE vcf_lines SET format = NULL where format = '.';
        UPDATE vcf_lines SET values = NULL where values = '.';
    """)

    db.sql("""
        ALTER TABLE
            vcf_lines
        ADD COLUMN
            vid VARCHAR;
    """)

    db.sql("""
        UPDATE
            vcf_lines
        SET
            vid = (
                chrom || ':' ||
                pos || '|' ||
                coalesce(ref, '.') || '>' ||
                coalesce(alt, '.')
            )
    """)

    db.sql("""
        INSERT INTO
            variants (
                vid,
                chrom,
                pos,
                id,
                ref,
                alt,
                qual,
                filters
            )
        SELECT
            vid,
            chrom,
            pos,
            id,
            ref,
            alt,
            qual,
            str_split(filters, ';')
        FROM
            vcf_lines;
    """)

    n_variants = db.table("vcf_lines").shape[0]
    max_batch_size = 100000
    n_batches = 1 + n_variants // max_batch_size
    batch_size = ceil(n_variants / n_batches)

    for i in range(n_batches):
        logging.info(f"Loading batch {i+1} of {n_batches}")
        offset = i * batch_size

        populate_vals(db, val_info_types, limit=batch_size, offset=offset)
        populate_info(db, val_info_types, limit=batch_size, offset=offset)

    db.sql("DROP TABLE IF EXISTS kv;")
    db.sql("DROP TABLE IF EXISTS kv_compound_info;")
    db.sql("DROP TABLE IF EXISTS vcf_lines;")
    db.unregister("sub_fields")

    for tbl in ["vals", "info"]:
        snake_case_col(db, tbl, "k")
        make_constraints(db, tbl)


def make_constraints(db: DuckDBPyConnection, tbl: str) -> None:
    logging.info(f"Applying constraints to {tbl}")

    db.sql(f"""
            INSERT INTO
                {tbl}_fk
            BY NAME
            SELECT
                *
            FROM
                {tbl};
        """)

    db.sql(f"DROP TABLE {tbl};")
    db.sql(f"ALTER TABLE {tbl}_fk RENAME TO {tbl};")


def snake_case_col(db: DuckDBPyConnection, tbl: str, col: str) -> None:
    logging.info(f"snake_casing {tbl}.{col}")

    snake_map = db.table(tbl)[col].distinct().df()
    snake_map[f"{col}_snake"] = snake_map[col].apply(snakecase)

    db.register("snake_map", snake_map)

    db.sql(f"""
        UPDATE
            {tbl}
        SET
            {col} = snake_map.{col}_snake
        FROM
            snake_map
        WHERE
            {tbl}.{col} = snake_map.{col};
    """)

    db.unregister("snake_map")


def populate_vals(
    db: DuckDBPyConnection,
    val_info_types: pd.DataFrame,
    limit: int = 0,
    offset: int = 0,
) -> None:
    db.sql("TRUNCATE kv;")

    db.sql(f"""
        INSERT INTO
            kv (
                vid,
                k,
                v
            )
        SELECT
            vid,
            unnest(str_split(format, ':')) as k,
            unnest(str_split(values, ':')) as v
        FROM (
            SELECT
                vid,
                format,
                values
            FROM
                vcf_lines
            ORDER BY
                chrom,
                pos,
                ref,
                alt
            LIMIT
                {limit}
            OFFSET
                {offset}
        );
    """)

    val_types = val_info_types.loc[
        val_info_types["kind"].eq("value")
        & val_info_types["id_db"].isin(db.table("kv")["k"].distinct().df()["k"])
    ].copy()

    cast_and_insert_v(
        db=db, src_table_name="kv", dest_table_name="vals", types_df=val_types
    )


def populate_info(
    db: DuckDBPyConnection,
    val_info_types: pd.DataFrame,
    limit: int = 0,
    offset: int = 0,
) -> None:
    db.sql("TRUNCATE kv;")

    info_types = val_info_types.loc[val_info_types["kind"].eq("info")].copy()
    info_types_expr = ", ".join(["'" + x + "'" for x in info_types["id_db"]])

    db.sql(f"""
        INSERT INTO
            kv (
                vid,
                k,
                v
            )
        SELECT
            vid,
            annot[1] as k,
            annot[2] as v
        FROM (
            SELECT
                vid,
                str_split(unnest(str_split(info, ';')), '=') AS annot
            FROM (
                SELECT
                    vid,
                    info
                FROM
                    vcf_lines
                ORDER BY
                    chrom,
                    pos,
                    ref,
                    alt
                LIMIT
                    {limit}
                OFFSET
                    {offset}
            )
        )
        WHERE
            k IN ({info_types_expr});
    """)

    info_types = info_types.loc[
        info_types["id_db"].isin(db.table("kv")["k"].distinct().df()["k"])
    ].copy()

    flag_fields = info_types["id_db"].loc[
        info_types["kind"].eq("info") & info_types["number"].eq("0")
    ]

    db.sql(f"""
        UPDATE
            kv
        SET
            v = 'true'
        WHERE
            k in ({', '.join(["'" + x + "'" for x in flag_fields])})
    """)

    simple_info_types = info_types.loc[~info_types["has_children"]].copy()
    compound_info_types = info_types.loc[info_types["has_children"]].copy()

    cast_and_insert_v(
        db=db, src_table_name="kv", dest_table_name="info", types_df=simple_info_types
    )

    compound_types_expr = ", ".join(
        ["'" + x + "'" for x in compound_info_types["id_db"]]
    )

    db.sql(f"""
        DELETE FROM
            kv
        WHERE
            k NOT IN ({compound_types_expr});
    """)

    db.sql(f"""
        INSERT INTO
            kv_compound_info (
                vid,
                k,
                k_sub,
                ix,
                v
            )
        SELECT
           vid,
           k,
           k_sub,
           ix,
           v
        FROM (
            WITH compound_info_split AS (
                SELECT
                    vid,
                    k,
                    str_split_regex(v, '\\s*\\|\\s*') AS v
                FROM
                    kv
            ),
            compound_info_exploded AS (
                SELECT
                    vid,
                    k,
                    unnest(v) AS v,
                    generate_subscripts(v, 1) AS ix
                FROM
                    compound_info_split
            )
            SELECT
                vid,
                compound_info_exploded.k,
                sub_fields.id_db as k_sub,
                sub_fields.ix,
                CASE
                    WHEN v IN ('', '.')
                THEN
                    NULL
                ELSE
                    v
                END AS v
            FROM
                compound_info_exploded
            INNER JOIN
                sub_fields
            ON
                compound_info_exploded.k = sub_fields.k
                AND
                compound_info_exploded.ix = sub_fields.ix
        );
    """)

    snake_case_col(db, tbl="kv_compound_info", col="k_sub")

    db.sql("TRUNCATE kv;")

    db.sql(f"""
        INSERT INTO
            kv (
                vid,
                k,
                v
            )
        SELECT
            vid,
            k,
            v
        FROM (
            WITH compound_info_maps AS (
                SELECT
                    vid,
                    k,
                    MAP(
                        list(k_sub ORDER BY ix),
                        list(v ORDER BY ix)
                    )::JSON::VARCHAR AS v
                FROM
                    kv_compound_info
                GROUP BY
                    vid,
                    k
            )
            SELECT
                vid,
                k,
                '[' || string_agg(v, ',') || ']' AS v
            FROM
                compound_info_maps
            GROUP BY
                vid,
                k
        )
    """)

    compound_info_types["is_list"] = False

    cast_and_insert_v(
        db=db,
        src_table_name="kv",
        dest_table_name="info",
        types_df=compound_info_types,
    )


def cast_and_insert_v(
    db: DuckDBPyConnection,
    src_table_name: str,
    dest_table_name: str,
    types_df: pd.DataFrame,
) -> None:
    for (v_col_name, is_list), g in types_df.groupby(["v_col_name", "is_list"]):
        k_ids_expr = ", ".join(["'" + x + "'" for x in g["id_db"]])

        if is_list:
            db.sql(f"""
                INSERT INTO
                    {dest_table_name} (
                        vid,
                        k,
                        {v_col_name}
                    )
                SELECT
                    vid,
                    k,
                    str_split(v, ',') as {v_col_name}
                FROM
                    {src_table_name}
                WHERE
                    k IN ({k_ids_expr});
            """)
        else:
            db.sql(f"""
                INSERT INTO
                    {dest_table_name} (
                        vid,
                        k,
                        {v_col_name}
                    )
                SELECT
                    vid,
                    k,
                    v
                FROM
                    {src_table_name}
                WHERE
                    k IN ({k_ids_expr});
            """)


def urldecode_cols(
    df: pd.DataFrame, url_encoded_col_name_regex: re.Pattern
) -> pd.DataFrame:
    col_has_percent = (
        df.select_dtypes("string").map(lambda x: "%" in str(x)).any(axis=0)
    )
    obs_percent_cols = col_has_percent[col_has_percent].index

    url_encoded_col_names = df.columns[df.columns.str.match(url_encoded_col_name_regex)]

    if not set(url_encoded_col_names).issuperset(set(obs_percent_cols)):
        # if this happens, we might need another CLI option to specify cols that have
        # percent signs but aren't actually URL-encoded
        others = set(obs_percent_cols).difference(set(url_encoded_col_names))
        logging.warning(f"Check VCF for additional URL-encoded info in {others}")

    df.loc[:, list(url_encoded_col_names)] = df.loc[:, list(url_encoded_col_names)].map(
        lambda x: unquote(x) if x is not pd.NA else pd.NA
    )

    return df
