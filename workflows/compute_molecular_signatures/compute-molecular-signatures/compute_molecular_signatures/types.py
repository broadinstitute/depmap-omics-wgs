from typing import TypeVar

import pandas as pd
import pandera.pandas as pa
from pandera.api.dataframe.model_config import BaseConfig
from pandera.typing import DataFrame, Series


class CoercedDataFrame(pa.DataFrameModel):
    class Config(BaseConfig):  # pyright: ignore
        coerce = True


class MutSigVariants(CoercedDataFrame):
    sample_id: Series[pd.StringDtype]
    chrom: Series[pd.StringDtype]
    pos: Series[pd.Int64Dtype]
    ref: Series[pd.StringDtype] = pa.Field(str_matches=r"[ACGT]+")
    alt: Series[pd.StringDtype] = pa.Field(str_matches=r"[ACGT]+")
    alt_count: Series[pd.Int64Dtype]
    af: Series[pd.Float64Dtype]
    vep_gnom_ade_af: Series[pd.Float64Dtype] = pa.Field(nullable=True)
    vep_gnom_adg_af: Series[pd.Float64Dtype] = pa.Field(nullable=True)


class Maf(CoercedDataFrame):
    Hugo_Symbol: Series[pd.StringDtype]
    Tumor_Sample_Barcode: Series[pd.StringDtype]
    Chromosome: Series[pd.StringDtype] = pa.Field(
        isin=[*[str(x) for x in range(1, 23)], "X", "Y"]
    )
    Start_Position: Series[pd.Int64Dtype]
    Reference_Allele: Series[pd.StringDtype] = pa.Field(str_matches=r"[ACGT\-]+")
    Tumor_Seq_Allele2: Series[pd.StringDtype] = pa.Field(str_matches=r"[ACGT\-]+")


PanderaBaseSchema = TypeVar("PanderaBaseSchema", bound=CoercedDataFrame)
TypedDataFrame = DataFrame
