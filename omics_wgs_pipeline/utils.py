from typing import Callable, Iterable, Type, TypeVar

import pandas as pd
from pydantic import BaseModel

from omics_wgs_pipeline.types import TypedDataFrame
from omics_wgs_pipeline.validators import CoercedDataFrame

B = TypeVar("B", bound=BaseModel)
T = TypeVar("T", bound=CoercedDataFrame)


def type_data_frame(
    df: pd.DataFrame, pandera_schema: Type[T], remove_unknown_cols: bool = True
) -> TypedDataFrame[T]:
    """
    Coerce a data frame into one specified by a Pandera schema and optionally remove
    unknown columns.

    :param df: a data frame
    :param pandera_schema: a Pandera schema
    :param remove_unknown_cols: remove columns not specified in the schema
    :return: a data frame validated with the provided Pandera schema
    """

    if remove_unknown_cols:
        df_cols = pandera_schema.to_schema().columns.keys()
        df = df.loc[:, df_cols].drop_duplicates()

    return TypedDataFrame[pandera_schema](df)


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
    `name_columns_with_parent` is `True`
    :return: a widened data frame
    """

    flattened_dict = {}

    for c, s in df.items():
        if isinstance(s.iloc[0], dict):
            # if the column contains dictionaries, recursively flatten them
            nested_df = pd.json_normalize(s.tolist())

            nested_df.columns = [
                sep.join([parent_key, c, str(col)])
                if name_columns_with_parent and parent_key != ""
                else str(col)
                for col in nested_df.columns
            ]

            flattened_dict.update(
                expand_dict_columns(
                    nested_df,
                    sep=sep,
                    name_columns_with_parent=name_columns_with_parent,
                    parent_key=str(c),
                )
            )

        else:
            # if not a dictionary, add the column as it is
            flattened_dict[c] = s

    df = pd.DataFrame(flattened_dict)

    # make sure there are no duplicate column names after expansion
    col_name_counts = df.columns.value_counts()

    if col_name_counts.gt(1).any():
        dup_names = set(col_name_counts[col_name_counts.gt(1)].index)
        raise NameError(
            f"Column names {dup_names} are duplicated. "
            "Try calling `expand_dict_columns` with `name_columns_with_parent=True`."
        )

    return df


def df_to_model(df: pd.DataFrame, pydantic_schema: Type[B]) -> list[B]:
    """
    Convert a Pandas data frame to a Pydantic model.

    :param df: a data frame
    :param pydantic_schema: the Pydantic schema to cast the data frame to
    :return: a Pydantic model
    """

    return [pydantic_schema(**x) for x in df.to_dict(orient="records")]


def model_to_df(
    model: B,
    pandera_schema: Type[T],
    records_key: str = "records",
    remove_unknown_cols: bool = True,
    mutator: Callable[[pd.DataFrame], pd.DataFrame] = lambda _: _,
) -> TypedDataFrame[T]:
    """
    Dump a Pydantic model and convert it to a data frame typed by a Pandera schema.

    :param model: a Pydandict model containing a list of objects keyed by `records_key`
    :param pandera_schema: the Pandera schema to cast the model to
    :param records_key: the key/method name in `model` containing the records
    :param remove_unknown_cols: remove columns not specified in the schema
    :param mutator: an optional function to call on the data frame before typing (e.g.
    to rename columns)
    """

    records = model.model_dump()[records_key]
    df = pd.DataFrame(records)
    df = mutator(df)
    return type_data_frame(df, pandera_schema, remove_unknown_cols)


def anti_join(
    x: pd.DataFrame, y: pd.DataFrame, on: str | Iterable[str]
) -> pd.DataFrame:
    """
    Anti join two data frames.

    :param x: a base data frame
    :param y: a data frame to use for filtering
    :param on: the columns to anti-join on
    :return: a data frame
    """

    if len(y) == 0:
        return x

    # make a data frame of just the join columns
    dummy = y.loc[:, on]

    # convert to data frame if `on` was a single column
    if isinstance(dummy, pd.Series):
        dummy = dummy.to_frame()

    dummy.loc[:, "dummy_col"] = 1  # indicator variable

    # attempt to join left data frame (`x`) to the dummy data frame
    merged = x.merge(dummy, on=on, how="left")

    # keep only the non-matches
    return merged.loc[merged["dummy_col"].isna(), x.columns.tolist()]
