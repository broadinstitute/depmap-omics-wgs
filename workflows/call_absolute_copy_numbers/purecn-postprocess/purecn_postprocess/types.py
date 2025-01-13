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


class Loh(CoercedDataFrame):
    size: Series[pd.UInt64Dtype]
    c: Series[pd.UInt64Dtype]
    m: Series[pd.UInt64Dtype] = pa.Field(nullable=True)
    type: Series[pd.StringDtype] = pa.Field(nullable=True)
