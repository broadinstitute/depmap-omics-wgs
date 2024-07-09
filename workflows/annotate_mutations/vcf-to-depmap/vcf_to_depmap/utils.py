import re
from pathlib import Path
from pprint import pp
from typing import Callable
from urllib.parse import unquote

import numpy as np
import pandas as pd
from caseconverter import snakecase


def get_vcf_info_and_format_dtypes(
    path: Path, compound_info_fields: set[str]
) -> pd.DataFrame:
    header_lines = []

    with open(path, "r") as f:
        while True:
            line = f.readline()

            if line.startswith("##FORMAT") or line.startswith("##INFO"):
                header_lines.append(line.strip())
            elif not line.startswith("#"):
                break

    arr = []

    for x in header_lines:
        kind = re.search(r"^##(\w+)", x).group(1).lower()
        interior = re.search(r"<(.+)>$", x).group(1)
        parts = re.findall(r'([A-Za-z0-9_]+)=(".*?"|[^,]+)', interior)

        d = {k: v.strip('"') for k, v in parts}
        d["kind"] = kind if kind == "info" else "value"
        d["has_subfields"] = d["ID"] in compound_info_fields
        d["ID"] = "__".join([d["kind"], snakecase(d["ID"])])

        if d["has_subfields"]:
            subfields = d["Description"].split(r"|")
            subfields[0] = subfields[0].split(":")[-1]

            d["subfields"] = [
                snakecase(re.search(r"[\s']*([^']+)[\s']*$", x).group(1))
                for x in subfields
            ]

        arr.append(d)

    df = pd.DataFrame(arr)
    df.columns = df.columns.str.lower()

    df["type"] = df["type"].replace(
        {
            "Integer": "Int64",
            "Float": "Float64",
            "String": "string",
            "Character": "string",
            "Flag": "boolean",
        }
    )

    return df


def read_vcf(path: Path) -> pd.DataFrame:
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=[
            "chromosome",
            "position",
            "id",
            "ref",
            "alt",
            "qual",
            "filter",
            "info",
            "format",
            "values",
        ],
        usecols=[
            "chromosome",
            "position",
            "ref",
            "alt",
            "filter",
            "info",
            "format",
            "values",
        ],
        comment="#",
        nrows=100,
    )

    df[["chromosome", "ref", "alt", "filter"]] = df[
        ["chromosome", "ref", "alt", "filter"]
    ].astype("string")
    df["position"] = df["position"].astype("int64")

    if df[["chromosome", "position", "ref", "alt"]].duplicated().any():
        raise ValueError("Duplicate variants detected in the VCF file.")

    return df


def process_vcf(
    df: pd.DataFrame,
    info_and_format_dtypes: pd.DataFrame,
    drop_cols: set[str],
    url_encoded_col_name_regex: str | None,
) -> pd.DataFrame:
    cols_to_drop = list(drop_cols)
    df.drop(columns=cols_to_drop, errors="ignore", inplace=True)

    df = expand_info_and_value_cols(df)
    df.drop(columns=cols_to_drop, errors="ignore", inplace=True)

    df = fix_info_and_value_dtypes(df, info_and_format_dtypes)
    df.drop(columns=cols_to_drop, errors="ignore", inplace=True)

    df = urldecode_cols(df, url_encoded_col_name_regex)

    return df


def expand_info_and_value_cols(df):
    df["info"] = df["info"].apply(parse_vcf_info)
    df["value"] = df.apply(
        lambda x: dict(zip(x["format"].split(":"), x["values"].split(":"))), axis=1
    )

    df = expand_dict_columns(df, col_name_formatter=snakecase)

    return df.drop(columns=["format", "values"])


def parse_vcf_info(info: str) -> dict[str, str]:
    parts = info.split(";")
    kv = [x.split("=") for x in parts]
    return dict(zip([x[0] for x in kv], [x[1] if len(x) == 2 else None for x in kv]))


def fix_info_and_value_dtypes(df, info_and_format_dtypes):
    obs_info_formats = info_and_format_dtypes.loc[
        info_and_format_dtypes["id"].isin(df.columns)
    ]

    for _, r in obs_info_formats.iterrows():
        if r["type"] == "boolean":
            df[r["id"]] = df[r["id"]].astype("boolean")
            df[r["id"]] = df[r["id"]].fillna(False)

        else:
            expanded = df[r["id"]].str.split(",", expand=True)

            if expanded.shape[1] == 1:
                new_col_names = [r["id"]]
                df[r["id"]] = df[r["id"]].astype(r["type"])
            else:
                new_col_names = [
                    "_".join([r["id"], str(x + 1)])
                    for x in list(range(0, expanded.shape[1]))
                ]

                df[new_col_names] = expanded.values
                df[new_col_names] = df[new_col_names].astype(r["type"])

                df = df.drop(columns=[r["id"]])

            if r["has_subfields"]:
                for c in new_col_names:
                    df[c] = df[c].str.strip("[]()").str.split("|")
                    df.loc[~df[c].isna(), c] = df.loc[~df[c].isna(), c].apply(
                        lambda x: dict(zip(r["subfields"], x))
                    )
                    df = expand_dict_columns(df)

                    df[df.columns[df.columns.str.startswith(c)]] = df[
                        df.columns[df.columns.str.startswith(c)]
                    ].astype(r["type"])

    assert np.dtype("O") not in list(df.dtypes)

    df = df[
        [
            "chromosome",
            "position",
            "ref",
            "alt",
            "filter",
            *df.columns[df.columns.str.startswith("value__")].sort_values(),
            *df.columns[df.columns.str.startswith("info__")].sort_values(),
        ]
    ]

    return df


def urldecode_cols(df, url_encoded_col_name_regex):
    col_has_percent = df.apply(lambda x: x.str.contains("%"), axis=1).any(axis=0)
    percent_cols = col_has_percent[col_has_percent].index

    url_encoded_col_names = df.columns[df.columns.str.match(url_encoded_col_name_regex)]

    if not set(url_encoded_col_names).issuperset(set(percent_cols)):
        raise ValueError("Check VCF for additional URL-encoded info")

    df.loc[:, url_encoded_col_names] = df.loc[:, url_encoded_col_names].map(
        lambda x: unquote(x) if x is not pd.NA else pd.NA
    )

    return df


def expand_dict_columns(
    df: pd.DataFrame,
    sep: str = "__",
    name_columns_with_parent: bool = True,
    parent_key: str = "",
    col_name_formatter: Callable[[str], str] = lambda _: _,
) -> pd.DataFrame:
    """
    Recursively expand columns in a data frame containing dictionaries into separate
    columns.

    :param df: a data frame
    :param sep: a separator character to use between `parent_key` and its column names
    :param name_columns_with_parent: whether to "namespace" nested column names using
    their parents' column names
    :param parent_key: the name of the parent column, applicable only if
    `name_columns_with_parent` is `True` (for recursion)
    :param col_name_formatter: an optional function to format resulting column names
    :return: a widened data frame
    """

    flattened_dict = {}

    for c, s in df.items():
        if isinstance(s.iloc[0], dict):
            # if the column contains dictionaries, recursively flatten them
            nested_df = pd.json_normalize(s.tolist())

            if name_columns_with_parent:
                # e.g. if current column `c` is "foo" and the nested data contains a
                # field "bar", the resulting column name is "foo__bar"

                nested_df.columns = [
                    sep.join(
                        [
                            parent_key,
                            col_name_formatter(str(c)),
                            col_name_formatter(str(col)),
                        ]
                    )
                    if parent_key != ""
                    else sep.join(
                        [col_name_formatter(str(c)), col_name_formatter(str(col))]
                    )
                    for col in nested_df.columns
                ]

            # recurse on the nested data
            flattened_dict.update(
                expand_dict_columns(
                    nested_df,
                    sep=sep,
                    name_columns_with_parent=name_columns_with_parent,
                    parent_key=col_name_formatter(str(c)),
                )
            )

        else:
            # if not a dictionary, add the column as is
            flattened_dict[c] = s

    df = pd.DataFrame(flattened_dict)

    if parent_key == "":
        # make sure there are no duplicate column names after all expansion is done
        col_name_counts = df.columns.value_counts()

        if col_name_counts.gt(1).any():
            dup_names = set(
                col_name_counts[col_name_counts.gt(1)].index,  # pyright: ignore
            )
            raise NameError(
                f"Column names {dup_names} are duplicated. Try calling "
                "`expand_dict_columns` with `name_columns_with_parent=True`."
            )

    return df
