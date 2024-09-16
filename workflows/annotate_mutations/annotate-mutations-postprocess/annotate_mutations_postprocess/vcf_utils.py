import concurrent.futures
import logging
import os
import re
from csv import QUOTE_NONE
from io import StringIO
from pathlib import Path
from urllib.parse import unquote

import bgzip
import duckdb
import numpy as np
import pandas as pd
from caseconverter import snakecase
from duckdb import DuckDBPyConnection
from vcf_info_merger.merge import get_header_lines

from annotate_mutations_postprocess.utils import expand_dict_columns


def vcf_to_wide(
    vcf_path: Path,
    compound_info_fields: set[str],
    info_subfield_types: dict[str, dict[str, str]],
    info_cols_ignored: set[str],
    url_encoded_col_name_regex: re.Pattern,
):
    info_and_format_dtypes = get_vcf_info_and_format_dtypes(
        vcf_path, compound_info_fields, info_cols_ignored, info_subfield_types
    )

    db_path = os.path.join(Path(__file__).parent.parent, "data", "variants.duckdb")
    try:
        os.remove(db_path)
    except OSError:
        pass

    with duckdb.connect(db_path) as db:
        set_up_db(db, info_and_format_dtypes)

        logging.info(f"Reading {vcf_path}")
        populate_db(db, vcf_path, info_and_format_dtypes=info_and_format_dtypes)
        pass


def set_up_db(db: DuckDBPyConnection, info_and_format_dtypes: pd.DataFrame) -> None:
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
        )
    """)

    val_cols = info_and_format_dtypes["col_def"].loc[
        info_and_format_dtypes["kind"].eq("value")
    ]

    db.sql(f"""
        CREATE TABLE IF NOT EXISTS vals (
            vid VARCHAR REFERENCES variants(vid),
            {', '.join(val_cols)}
        )
    """)

    nested_info_fields = info_and_format_dtypes.loc[
        info_and_format_dtypes["kind"].eq("info")
        & info_and_format_dtypes["has_children"]
    ]

    for _, field in nested_info_fields.iterrows():
        sub_fields = info_and_format_dtypes["col_def"].loc[
            info_and_format_dtypes["parent_id_db"].eq(field["id_db"])
        ]

        db.sql(f"""
            CREATE TYPE info_{field['id_db']} AS STRUCT (
                {', '.join(sub_fields)}
            )
        """)

    info_cols = info_and_format_dtypes["col_def"].loc[
        info_and_format_dtypes["kind"].eq("info")
    ]

    db.sql(f"""
        CREATE TABLE IF NOT EXISTS info (
            vid VARCHAR REFERENCES variants(vid),
            {', '.join(info_cols)}
        )
    """)


def get_vcf_info_and_format_dtypes(
    path: Path,
    compound_info_fields: set[str],
    info_cols_ignored: set[str],
    info_subfield_types: dict[str, dict[str, str]],
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

        d["id_db"] = snakecase(d["id"])

        if d["id"] in compound_info_fields:
            d["has_children"] = True
            d["type"] = f"info_{d['id_db']}"

            desc = re.search(r"^.+:['\s]*([^']+)['\s]*$", d["description"]).group(1)
            subfields = re.split(r"\s*\|\s*", desc)

            for ix, s in enumerate(subfields):
                try:
                    sub_type = info_subfield_types[d["id"]][s]
                except KeyError:
                    sub_type = "VARCHAR"

                dsub = {
                    "id": s,
                    "id_db": snakecase(s),
                    "has_children": False,
                    "number": "1",
                    "type": sub_type,
                    "kind": "sub_info",
                    "parent_id": d["id"],
                    "parent_id_db": d["id_db"],
                    "ix": ix,
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

    df["col_def"] = df["id_db"] + " " + df["type"]
    df.loc[~df["number"].isin({"0", "1"}), "col_def"] += "[]"

    return df


def populate_db(
    db: DuckDBPyConnection,
    path: Path | bgzip.BGZipReader,
    info_and_format_dtypes: pd.DataFrame,
) -> None:
    memoryview_batch_size = 2**24
    chunk_size = 2**18

    with open(path, "rb") as raw:
        seen_header = False
        batch = ""

        with bgzip.BGZipReader(
            raw,
            buffer_size=memoryview_batch_size * 3,
            num_threads=1,
            raw_read_chunk_size=chunk_size,
        ) as f:
            with concurrent.futures.ThreadPoolExecutor() as executor:
                futures = []

                while True:
                    # read a block of bytes
                    if not (d := f.read(memoryview_batch_size)):
                        break

                    # concat the latest chunk of text
                    text = d.tobytes().decode()
                    d.release()
                    batch += text

                    if not seen_header:
                        m = re.search(r"\n(#CHROM[^\n]+)\n", batch, re.MULTILINE)

                        if m is not None:
                            logging.info("Found column header. Reading variants.")

                            # done with VCF header, confirm col header format
                            assert m[1].split("\t")[:-1] == [
                                "#CHROM",
                                "POS",
                                "ID",
                                "REF",
                                "ALT",
                                "QUAL",
                                "FILTER",
                                "INFO",
                                "FORMAT",
                            ]

                            batch = batch[m.end() :]
                            seen_header = True

                        else:
                            continue  # still looking for col header

                    batch_lines, batch = batch.rsplit("\n", 1)

                    logging.info(
                        "\t".join(batch_lines.split("\n", 1)[0].split("\t", 5)[:5])
                    )

                    futures.append(
                        executor.submit(
                            insert_batch, db, batch_lines, info_and_format_dtypes
                        )
                    )

                concurrent.futures.wait(futures)


def insert_batch(
    db: DuckDBPyConnection, batch_lines: str, info_and_format_dtypes: pd.DataFrame
) -> None:
    cur = db.cursor()

    df = pd.read_csv(
        StringIO(batch_lines),
        sep="\t",
        header=None,
        names=[
            "chrom",
            "pos",
            "id",
            "ref",
            "alt",
            "qual",
            "filters",
            "info",
            "format",
            "values",
        ],
        dtype="string",
        na_values=".",
        keep_default_na=False,
        quoting=QUOTE_NONE,
        encoding_errors="backslashreplace",
    )

    df["vid"] = df["chrom"] + ":" + df["pos"] + "|" + df["ref"] + ">" + df["alt"]

    df["pos"] = df["pos"].astype("int32")

    variants_df = df.loc[
        :, ["vid", "chrom", "pos", "id", "ref", "alt", "qual", "filters"]
    ]

    variants_df["filters"] = variants_df["filters"].apply(lambda x: x.split(";"))

    cur.sql("INSERT INTO variants BY NAME SELECT * FROM variants_df")

    vals_df = df.loc[:, ["vid", "format", "values"]].set_index("vid")

    vals_df = vals_df.apply(
        lambda x: dict(zip(x["format"].split(":"), x["values"].split(":"))), axis=1
    ).to_frame()

    vals_df = expand_dict_columns(
        vals_df, col_name_formatter=snakecase, name_columns_with_parent=False
    )

    multi_val_cols = info_and_format_dtypes["id_db"].loc[
        info_and_format_dtypes["id_db"].isin(vals_df.columns)
        & info_and_format_dtypes["kind"].eq("value")
        & ~info_and_format_dtypes["number"].isin({"0", "1"})
    ]

    for c in multi_val_cols:
        vals_df.loc[vals_df[c].notna(), c] = vals_df.loc[vals_df[c].notna(), c].apply(
            lambda x: x.split(",")
        )

    vals_df = vals_df.reset_index()

    cur.sql("INSERT INTO vals BY NAME SELECT * FROM vals_df")

    info_df = df.loc[:, ["vid", "info"]].set_index("vid")
    info_df["info"] = info_df["info"].apply(parse_vcf_info)

    info_df = expand_dict_columns(
        info_df, col_name_formatter=snakecase, name_columns_with_parent=False
    )

    info_cols = info_and_format_dtypes.loc[
        info_and_format_dtypes["id_db"].isin(info_df.columns)
        & info_and_format_dtypes["kind"].eq("info")
    ]

    info_df = info_df.loc[:, info_df.columns.isin(info_cols["id_db"])]

    multi_info_cols = info_cols["id_db"].loc[~info_cols["number"].isin({"0", "1"})]

    for c in multi_info_cols:
        info_df.loc[info_df[c].notna(), c] = info_df.loc[info_df[c].notna(), c].apply(
            lambda x: x.split(",")
        )

    nested_fields = info_cols.loc[info_cols["has_children"]]

    for _, info_field in nested_fields.iterrows():
        sub_col_names = info_and_format_dtypes.loc[
            info_and_format_dtypes["parent_id_db"].eq(info_field["id_db"])
        ].sort_values("ix")["id_db"]

        sub_info = (
            info_df[info_field["id_db"]]
            .dropna()
            .explode()
            .str.split(r"\s*\|\s*", expand=True)
            .replace({"": pd.NA, ".": pd.NA})
        )

        sub_info.columns = sub_col_names

        sub_info = pd.Series(
            sub_info.to_dict(orient="records"),
            index=sub_info.index,
            name=info_field["id_db"],
        )

        sub_info = sub_info.groupby(sub_info.index).agg(list)

        info_df.update(sub_info)
        info_df = info_df.rename(columns={info_field["id_db"]: info_field["id_db"]})

    info_df = info_df.reset_index()

    cur.sql("INSERT INTO info BY NAME SELECT * FROM info_df")


def clean_vcf(
    df: pd.DataFrame,
    info_and_format_dtypes: pd.DataFrame,
    bool_cols: set[str],
    drop_cols: set[str],
    url_encoded_col_name_regex: re.Pattern,
) -> pd.DataFrame:
    cols_to_drop = list(drop_cols)
    df.drop(columns=cols_to_drop, errors="ignore", inplace=True)

    df = expand_and_cast(df, info_and_format_dtypes)
    df.drop(columns=cols_to_drop, errors="ignore", inplace=True)

    df.replace({"": pd.NA, ".": pd.NA}, inplace=True)

    logging.info("URL-decoding columns")
    df = urldecode_cols(df, url_encoded_col_name_regex)

    logging.info("Casting boolean columns")
    df = convert_booleans(df, bool_cols)

    return df


def parse_vcf_info(info: str) -> dict[str, str | None]:
    parts = info.split(";")
    kv = [x.split("=") for x in parts]
    return dict(zip([x[0] for x in kv], [x[1] if len(x) == 2 else None for x in kv]))


def expand_and_cast(
    info_df: pd.DataFrame, info_and_format_dtypes: pd.DataFrame
) -> pd.DataFrame:
    obs_info_formats = info_and_format_dtypes.loc[
        info_and_format_dtypes["kind"].eq("info")
        & info_and_format_dtypes["id"].isin(info_df.columns)
    ]

    for _, r in obs_info_formats.iterrows():
        if r["type"] == "boolean":
            info_df[r["id"]] = info_df[r["id"]].astype("boolean")
            info_df[r["id"]] = info_df[r["id"]].fillna(False)

        else:
            expanded = info_df[r["id"]].str.split(",", expand=True)

            if expanded.shape[1] == 1:
                new_col_names = [r["id"]]
            else:
                new_col_names = [
                    r["id"],
                    *[
                        ".".join([r["id"], str(x + 1)])
                        for x in list(range(1, expanded.shape[1]))
                    ],
                ]

                info_df[new_col_names] = expanded.values

            if r["has_subfields"]:
                logging.info(f"Expanding {r['id']} subfields")

                for c in new_col_names:
                    if c == "info__ann":
                        print(0)
                    info_df[c] = info_df[c].str.strip("[]()").str.split("|")
                    info_df.loc[~info_df[c].isna(), c] = info_df.loc[
                        ~info_df[c].isna(), c
                    ].apply(lambda x: dict(zip(r["subfields"], x)))
                    info_df = expand_dict_columns(info_df)

                    regex = re.compile(f"^{re.escape(c)}__")
                    new_col_names_w_subfields = info_df.columns[
                        info_df.columns.str.match(regex)
                    ]

                    info_df[new_col_names_w_subfields] = info_df[
                        new_col_names_w_subfields
                    ].astype(r["type"])

                    if len(new_col_names) > 1:
                        digit_move_map = dict(
                            zip(
                                new_col_names_w_subfields,
                                (
                                    r["id"]
                                    + "__"
                                    + new_col_names_w_subfields.str.replace(
                                        regex, "", regex=True
                                    )
                                    + "."
                                    + re.search(r"\d+$", c)[0]
                                ),
                            )
                        )

                        info_df = info_df.rename(columns=digit_move_map)

            else:
                info_df[new_col_names] = info_df[new_col_names].astype(r["type"])

    assert np.dtype("O") not in list(info_df.dtypes)

    info_df = info_df[
        [
            "chromosome",
            "position",
            "ref",
            "alt",
            *info_df.columns[info_df.columns.str.startswith("value__")].sort_values(),
            *info_df.columns[info_df.columns.str.startswith("filter__")].sort_values(),
            *info_df.columns[info_df.columns.str.startswith("info__")].sort_values(),
        ]
    ]

    return info_df


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


def convert_booleans(df: pd.DataFrame, bool_cols: set[str]) -> pd.DataFrame:
    df_strings = df.select_dtypes(include="string")
    df_strings = df_strings.loc[:, df_strings.notna().any(axis=0)]

    true_vals = {"True", "true", "Y", "y"}
    false_vals = {"False", "false", "N", "n"}

    col_is_boollike = df_strings.isin({*true_vals, *false_vals, pd.NA}).all(axis=0)

    obs_bool_cols = col_is_boollike[col_is_boollike].index

    if not set(obs_bool_cols).issubset(bool_cols):
        others = set(obs_bool_cols).difference(bool_cols)
        logging.warning(f"Check VCF for additional booleans: {others}")

    df[obs_bool_cols] = df[obs_bool_cols].astype("object")

    for c in obs_bool_cols:
        df.loc[df[c].isin({"True", "true", "Y", "y"}), c] = True
        df.loc[df[c].isin({"False", "false", "N", "n"}), c] = False

    df[obs_bool_cols] = df[obs_bool_cols].astype("boolean")

    return df
