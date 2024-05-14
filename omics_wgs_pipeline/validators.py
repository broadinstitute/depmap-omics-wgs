import pandas as pd
import pandera as pa
from pandera.api.pandas.model_config import BaseConfig as PaBaseConfig
from pandera.typing import Series


class CoercedDataFrame(pa.DataFrameModel):
    class Config(PaBaseConfig):
        coerce = True  # convert to indicated dtype upon TypedDataFrame init


class SampleToPreprocess(CoercedDataFrame):
    sample_id: Series[pd.StringDtype]
    delivery_file_format: Series[pd.StringDtype] = pa.Field(isin={"CRAM", "BAM"})
    delivery_cram_bam: Series[pd.StringDtype]
    delivery_crai_bai: Series[pd.StringDtype]
    delivery_ref_fasta: Series[pd.StringDtype]
    delivery_ref_fasta_index: Series[pd.StringDtype]
