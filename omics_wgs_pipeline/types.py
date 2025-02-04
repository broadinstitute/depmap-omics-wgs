from typing import Any, Optional, TypeVar

import httpx
import pandas as pd
import pandera as pa
import pandera.typing
from nebelung.types import CoercedDataFrame
from pandera.typing import Series
from pydantic import BaseModel

from gumbo_gql_client.gumbo_client import GumboClient as AriadneGumboClient


class GumboClient(AriadneGumboClient):
    def __init__(
        self,
        url: str,
        username: str,
        headers: dict[str, str],
        http_client: Optional[httpx.Client] = None,
    ):
        super().__init__(url=url, headers=headers, http_client=http_client)
        self.username = username  # store username on this object for use in mutations


class GumboWgsSequencing(CoercedDataFrame):
    omics_sequencing_id: Series[pd.StringDtype]
    sequencing_alignment_source: Series[pd.StringDtype] = pa.Field(isin=["GP", "CDS"])
    reference_genome: Series[pd.StringDtype]
    url: Series[pd.StringDtype] = pa.Field(unique=True)
    index_url: Series[pd.StringDtype] = pa.Field(unique=True)


class TerraSample(CoercedDataFrame):
    sample_id: Series[pd.StringDtype] = pa.Field(unique=True)
    delivery_cram_bam: Series[pd.StringDtype] = pa.Field(nullable=True)
    delivery_crai_bai: Series[pd.StringDtype] = pa.Field(nullable=True)
    delivery_file_format: Series[pd.StringDtype] = pa.Field(
        isin={"CRAM", "BAM"}, nullable=True
    )
    delivery_ref: Series[pd.StringDtype]
    delivery_ref_amb: Series[pd.StringDtype]
    delivery_ref_ann: Series[pd.StringDtype]
    delivery_ref_bwt: Series[pd.StringDtype]
    delivery_ref_dict: Series[pd.StringDtype]
    delivery_ref_fasta: Series[pd.StringDtype]
    delivery_ref_fasta_index: Series[pd.StringDtype]
    delivery_ref_pac: Series[pd.StringDtype]
    delivery_ref_sa: Series[pd.StringDtype]
    analysis_ready_bam: Series[pd.StringDtype] = pa.Field(nullable=True)
    analysis_ready_bai: Series[pd.StringDtype] = pa.Field(nullable=True)
    ref: Series[pd.StringDtype]
    ref_amb: Series[pd.StringDtype]
    ref_ann: Series[pd.StringDtype]
    ref_bwt: Series[pd.StringDtype]
    ref_dict: Series[pd.StringDtype]
    ref_fasta: Series[pd.StringDtype]
    ref_fasta_index: Series[pd.StringDtype]
    ref_pac: Series[pd.StringDtype]
    ref_sa: Series[pd.StringDtype]


class GumboTaskEntity(CoercedDataFrame):
    id: Series[pd.StringDtype]
    sequencing_id: Series[pd.StringDtype] = pa.Field(nullable=True)


class GumboTaskResult(CoercedDataFrame):
    sample_id: Series[pd.StringDtype]
    label: Series[pd.StringDtype]
    url: Series[pd.StringDtype] = pa.Field(nullable=True)
    value: Series[dict[str, Any]] = pa.Field(nullable=True)
    workflow_name: Series[pd.StringDtype] = pa.Field(nullable=True)
    workflow_version: Series[pd.StringDtype] = pa.Field(nullable=True)
    completed_at: Series[pd.Timestamp]


class GcsObject(CoercedDataFrame):
    url: Series[pd.StringDtype]
    size: Series[pd.Int64Dtype]
    crc32c_hash: Series[pd.StringDtype]
    created_at: Series[pd.Timestamp]


PydanticBaseModel = TypeVar("PydanticBaseModel", bound=BaseModel)
PanderaBaseSchema = TypeVar("PanderaBaseSchema", bound=CoercedDataFrame)
TypedDataFrame = pandera.typing.DataFrame
T = TypeVar("T")
