import re
from pathlib import Path
from pprint import pp

import numpy as np
import pandas as pd


def get_vcf_info_and_format_dtypes(path: Path) -> pd.DataFrame:
    header_lines = []

    with open(path, "r") as f:
        while True:
            line = f.readline()

            if line.startswith("##FORMAT") or line.startswith("##INFO"):
                header_lines.append(line.strip())
            elif not line.startswith("#"):
                break

    rows = []

    for x in header_lines:
        kind = re.search(r"^##(\w+)", x).group(1).lower()
        interior = re.search(r"<([^>]+)>$", x).group(1)
        parts = re.findall(r'([A-Z a-z 0-9 _]+)=(".*?"|[^,]+)', interior)
        kv = {k: v.strip('"') for k, v in parts}
        kv["kind"] = kind
        rows.append(kv)

    df = pd.DataFrame(rows)
    df = df.drop_duplicates()
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


def read_vcf(path: Path, info_and_format_dtypes: pd.DataFrame) -> pd.DataFrame:
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

    df = expand_info_and_value_cols(df)
    obs_info_formats = collect_castable_info_and_values(df, info_and_format_dtypes)
    df = fix_info_and_value_dtypes(df, obs_info_formats)

    return df


def expand_info_and_value_cols(df):
    df["info"] = df["info"].apply(parse_vcf_info)
    df["value"] = df.apply(
        lambda x: dict(zip(x["format"].split(":"), x["values"].split(":"))), axis=1
    )

    df = expand_dict_columns(df)

    return df.drop(columns=["format", "values"])


def parse_vcf_info(info: str) -> dict[str, str]:
    parts = info.split(";")
    kv = [x.split("=") for x in parts]
    return dict(zip([x[0] for x in kv], [x[1] if len(x) == 2 else None for x in kv]))


def collect_castable_info_and_values(df, info_and_format_dtypes):
    format_rows = info_and_format_dtypes.loc[
        info_and_format_dtypes["kind"].eq("format")
    ].copy()

    info_rows = info_and_format_dtypes.loc[
        info_and_format_dtypes["kind"].eq("info")
    ].copy()

    format_rows["id"] = "value__" + format_rows["id"]
    info_rows["id"] = "info__" + info_rows["id"]

    info_and_formats = pd.concat([info_rows, format_rows], ignore_index=True)
    return info_and_formats.loc[info_and_formats["id"].isin(df.columns)]


def fix_info_and_value_dtypes(df, obs_info_formats):
    for _, r in obs_info_formats.iterrows():
        if r["type"] == "boolean":
            df[r["id"]] = df[r["id"]].astype("boolean")
            df[r["id"]] = df[r["id"]].fillna(False)

        else:
            expanded = df[r["id"]].str.split(",", expand=True)

            if expanded.shape[1] == 1:
                df[r["id"]] = df[r["id"]].astype(r["type"])
            else:
                new_col_names = [
                    "_".join([r["id"], str(x + 1)])
                    for x in list(range(0, expanded.shape[1]))
                ]

                df[new_col_names] = expanded.values
                df[new_col_names] = df[new_col_names].astype(r["type"])

                df = df.drop(columns=[r["id"]])

    assert np.dtype("O") not in list(df.dtypes)

    return df


def expand_dict_columns(
    df: pd.DataFrame,
    sep: str = "__",
    name_columns_with_parent: bool = True,
    parent_key: str = "",
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
                    sep.join([parent_key, str(c), str(col)])
                    if parent_key != ""
                    else sep.join([str(c), str(col)])
                    for col in nested_df.columns
                ]

            # recurse on the nested data
            flattened_dict.update(
                expand_dict_columns(
                    nested_df,
                    sep=sep,
                    name_columns_with_parent=name_columns_with_parent,
                    parent_key=str(c),
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
