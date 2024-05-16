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
    bam: Series[pd.StringDtype] = pa.Field(nullable=True)
    bai: Series[pd.StringDtype] = pa.Field(nullable=True)


class GumboWgsSequencing(CoercedDataFrame):
    sequencing_id: Series[pd.StringDtype]
    hg19_bai_filepath: Series[pd.StringDtype] = pa.Field(nullable=True)
    hg19_bam_filepath: Series[pd.StringDtype] = pa.Field(nullable=True)
    hg38_crai_filepath: Series[pd.StringDtype] = pa.Field(nullable=True)
    hg38_cram_filepath: Series[pd.StringDtype] = pa.Field(nullable=True)
    bai_filepath: Series[pd.StringDtype] = pa.Field(nullable=True)
    bam_filepath: Series[pd.StringDtype] = pa.Field(nullable=True)


class GumboTaskResult(CoercedDataFrame):
    sample_id: Series[pd.StringDtype]
    label: Series[pd.StringDtype]
    url: Series[pd.StringDtype] = pa.Field(nullable=True)
