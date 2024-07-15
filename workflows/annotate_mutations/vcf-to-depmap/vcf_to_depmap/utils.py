import re
from typing import Callable, Iterable

import pandas as pd


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
        fvi = s.first_valid_index()

        if fvi is not None and isinstance(s.loc[fvi], dict):
            # if the column contains dictionaries, recursively flatten them
            nested_df = pd.json_normalize(s.tolist())
            nested_df.index = df.index

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


def cs(df: pd.DataFrame, prefix: str) -> list[str]:
    regex = re.compile(rf"^{re.escape(prefix)}(_\d+)?$")
    return list(df.columns[df.columns.str.match(regex)])
