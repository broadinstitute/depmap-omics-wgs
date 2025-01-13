from typing import TypeVar

import pandas as pd
import pandera as pa
import pandera.typing
from pandera.api.pandas.model_config import BaseConfig as PaBaseConfig
from pandera.typing import Series


class CoercedDataFrame(pa.DataFrameModel):
    class Config(PaBaseConfig):
        coerce = True  # convert to indicated dtype upon TypedDataFrame init


PanderaBaseSchema = TypeVar("PanderaBaseSchema", bound=CoercedDataFrame)
TypedDataFrame = pandera.typing.DataFrame


class Solution(CoercedDataFrame):
    purity: Series[pd.Float64Dtype]
    ploidy: Series[pd.Float64Dtype]
    contamination: Series[pd.Float64Dtype]
    flagged: Series[pd.BooleanDtype]
    curated: Series[pd.BooleanDtype]
    comment: Series[pd.StringDtype] = pa.Field(nullable=True)


class Loh(CoercedDataFrame):
    size: Series[pd.UInt64Dtype]
    c: Series[pd.UInt64Dtype]
    m: Series[pd.UInt64Dtype] = pa.Field(nullable=True)
    type: Series[pd.StringDtype] = pa.Field(nullable=True)
