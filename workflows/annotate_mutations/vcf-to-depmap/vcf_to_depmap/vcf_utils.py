import re
from pathlib import Path
from urllib.parse import unquote

import numpy as np
import pandas as pd
from _csv import QUOTE_NONE
from caseconverter import snakecase

from vcf_to_depmap.utils import expand_dict_columns


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
        dtype="string",
        na_values=".",
        keep_default_na=False,
        quoting=QUOTE_NONE,
        low_memory=False,
        # nrows=1000,
        encoding_errors="backslashreplace",
    )

    header_line = (
        df.loc[df["chromosome"].eq("#CHROM")]
        .drop(columns="values")
        .to_dict(orient="records")
    )[0]

    assert header_line == {
        "chromosome": "#CHROM",
        "position": "POS",
        "id": "ID",
        "ref": "REF",
        "alt": "ALT",
        "qual": "QUAL",
        "filter": "FILTER",
        "info": "INFO",
        "format": "FORMAT",
    }

    df = df.loc[~df["chromosome"].str.startswith("#")]
    df = df.drop(columns=["id", "qual"])
    df["position"] = df["position"].astype("int64")

    assert ~df[["chromosome", "position", "ref", "alt"]].duplicated().any()
    assert ~df.isna().any().any()

    return df.reset_index(drop=True)


def clean_vcf(
    df: pd.DataFrame,
    info_and_format_dtypes: pd.DataFrame,
    drop_cols: set[str],
    na_cols: set[str],
    bool_cols: set[str],
    url_encoded_col_name_regex: re.Pattern,
    funco_sanitized_col_name_regex: re.Pattern,
) -> pd.DataFrame:
    cols_to_drop = list(drop_cols)
    df.drop(columns=cols_to_drop, errors="ignore", inplace=True)

    df = expand_info_value_filters(df)
    df.drop(columns=cols_to_drop, errors="ignore", inplace=True)

    df = expand_and_cast(df, info_and_format_dtypes)
    df.drop(columns=cols_to_drop, errors="ignore", inplace=True)

    df.replace({"": pd.NA}, inplace=True)
    df.replace({"": pd.NA, ".": pd.NA}, inplace=True)

    df = urldecode_cols(df, url_encoded_col_name_regex, funco_sanitized_col_name_regex)
    # df = remove_na_cols(df, na_cols)
    df = convert_booleans(df, bool_cols)

    return df


def expand_info_value_filters(df: pd.DataFrame) -> pd.DataFrame:
    df["info"] = df["info"].apply(parse_vcf_info)
    df["value"] = df.apply(
        lambda x: dict(zip(x["format"].split(":"), x["values"].split(":"))), axis=1
    )
    df = expand_filters(df)
    df = expand_dict_columns(df, col_name_formatter=snakecase)

    return df.drop(columns=["format", "values"])


def parse_vcf_info(info: str) -> dict[str, str]:
    parts = info.split(";")
    kv = [x.split("=") for x in parts]
    return dict(zip([x[0] for x in kv], [x[1] if len(x) == 2 else None for x in kv]))


def expand_filters(df: pd.DataFrame) -> pd.DataFrame:
    filters_long = df["filter"].str.split(";").explode().to_frame()
    obs_filters = list(filters_long["filter"].unique())

    filters_long[obs_filters] = False

    for f in obs_filters:
        filters_long[f] = filters_long["filter"].eq(f)

    filters_long = filters_long.drop(columns="filter")
    filters_long.columns = "filter__" + filters_long.columns.str.lower()
    filters_long = filters_long.groupby(filters_long.index)[filters_long.columns].any()
    filters_long = filters_long.astype("boolean")

    return df.join(filters_long, how="left").drop(columns="filter")


def expand_and_cast(df, info_and_format_dtypes):
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
            else:
                new_col_names = [
                    "_".join([r["id"], str(x + 1)])
                    for x in list(range(0, expanded.shape[1]))
                ]

                df[new_col_names] = expanded.values
                df = df.drop(columns=[r["id"]])

            if r["has_subfields"]:
                for c in new_col_names:
                    df[c] = df[c].str.strip("[]()").str.split("|")
                    df.loc[~df[c].isna(), c] = df.loc[~df[c].isna(), c].apply(
                        lambda x: dict(zip(r["subfields"], x))
                    )
                    df = expand_dict_columns(df)

                    new_col_names_w_subfields = df.columns[df.columns.str.startswith(c)]

                    df[new_col_names_w_subfields] = df[
                        new_col_names_w_subfields
                    ].astype(r["type"])

                    if len(new_col_names) > 1:
                        regex = re.compile(f"^{c}__")

                        digit_move_map = dict(
                            zip(
                                new_col_names_w_subfields,
                                (
                                    r["id"]
                                    + "__"
                                    + new_col_names_w_subfields.str.replace(
                                        regex, "", regex=True
                                    )
                                    + "_"
                                    + re.search(r"\d+$", c)[0]
                                ),
                            )
                        )

                        df = df.rename(columns=digit_move_map)

            else:
                df[new_col_names] = df[new_col_names].astype(r["type"])

    assert np.dtype("O") not in list(df.dtypes)

    df = df[
        [
            "chromosome",
            "position",
            "ref",
            "alt",
            *df.columns[df.columns.str.startswith("value__")].sort_values(),
            *df.columns[df.columns.str.startswith("filter__")].sort_values(),
            *df.columns[df.columns.str.startswith("info__")].sort_values(),
        ]
    ]

    return df


def urldecode_cols(
    df: pd.DataFrame,
    url_encoded_col_name_regex: re.Pattern,
    funco_sanitized_col_name_regex: re.Pattern,
) -> pd.DataFrame:
    col_has_percent = df.apply(lambda x: x.str.contains("%"), axis=1).any(axis=0)
    obs_percent_cols = col_has_percent[col_has_percent].index

    url_encoded_col_names = df.columns[df.columns.str.match(url_encoded_col_name_regex)]
    funco_sanitized_col_names = df.columns[
        df.columns.str.match(funco_sanitized_col_name_regex)
    ]
    either_col_names = set(url_encoded_col_names).union(set(funco_sanitized_col_names))

    if not set(either_col_names).issuperset(set(obs_percent_cols)):
        # if this happens, we might need another CLI option to specify cols that have
        # percent signs but aren't actually URL-encoded
        others = set(obs_percent_cols).difference(set(either_col_names))
        raise ValueError(f"Check VCF for additional URL-encoded info in {others}")

    for c in funco_sanitized_col_names:
        df[c] = df[c].str.replace(
            r"_(%[A-Za-z0-9]{2})_", lambda x: x.group(1), regex=True
        )

    df.loc[:, list(either_col_names)] = df.loc[:, list(either_col_names)].map(
        lambda x: unquote(x) if x is not pd.NA else pd.NA
    )

    return df


def remove_na_cols(df: pd.DataFrame, na_cols: set[str]) -> pd.DataFrame:
    col_is_na = df.isna().all(axis=0)
    obs_na_cols = col_is_na[col_is_na].index
    # assert set(obs_na_cols).issubset(na_cols)

    return df.drop(columns=col_is_na[col_is_na].index)


def convert_booleans(df: pd.DataFrame, bool_cols: set[str]) -> pd.DataFrame:
    df_strings = df.select_dtypes(include="string")

    true_vals = {"True", "true", "Y", "y"}
    false_vals = {"False", "false", "N", "n"}

    col_is_boollike = df_strings.isin({*true_vals, *false_vals, pd.NA}).all(axis=0)

    obs_bool_cols = col_is_boollike[col_is_boollike].index
    # assert set(obs_bool_cols) == bool_cols

    df[obs_bool_cols] = df[obs_bool_cols].astype("object")

    for c in obs_bool_cols:
        df.loc[df[c].isin({"True", "true", "Y", "y"}), c] = True
        df.loc[df[c].isin({"False", "false", "N", "n"}), c] = False

    df[obs_bool_cols] = df[obs_bool_cols].astype("boolean")

    return df
