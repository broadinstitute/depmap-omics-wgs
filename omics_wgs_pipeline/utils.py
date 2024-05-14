from typing import Type, TypeVar

import pandas as pd

from omics_wgs_pipeline.types import TypedDataFrame
from omics_wgs_pipeline.validators import CoercedDataFrame

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
    df: pd.DataFrame, parent_key: str = "", sep: str = "_"
) -> pd.DataFrame:
    """
    Recursively expand columns in a data frame containing dictionaries into separate
    columns.

    :param df: a data frame
    :param parent_key: a prefix to use for nested column names
    :param sep: a separator character to use between `parent_key` and its column names

    :return: a widened data frame
    """

    flattened_dict = {}

    for c, s in df.items():
        if isinstance(s.iloc[0], dict):
            # if the column contains dictionaries, recursively flatten them
            nested_df = pd.json_normalize(s.tolist())

            nested_df.columns = [
                f"{parent_key}{sep}{col}" if parent_key else col
                for col in nested_df.columns
            ]

            flattened_dict.update(
                expand_dict_columns(nested_df, parent_key=str(c), sep=sep)
            )

        else:
            # if not a dictionary, add the column as it is
            flattened_dict[c] = s

    return pd.DataFrame(flattened_dict)
