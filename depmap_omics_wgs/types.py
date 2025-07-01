from typing import Optional, TypeVar

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
    model_id: Series[pd.StringDtype]
    cell_line_name: Series[pd.StringDtype]
    stripped_cell_line_name: Series[pd.StringDtype]
    model_condition_id: Series[pd.StringDtype]
    omics_profile_id: Series[pd.StringDtype]
    omics_sequencing_id: Series[pd.StringDtype]
    sequencing_alignment_source: Series[pd.StringDtype] = pa.Field(isin=["GP", "CDS"])
    reference_genome: Series[pd.StringDtype]
    url: Series[pd.StringDtype] = pa.Field(unique=True)
    index_url: Series[pd.StringDtype] = pa.Field(unique=True)
    size: Series[pd.Int64Dtype] = pa.Field(unique=True)


class TerraSample(CoercedDataFrame):
    sample_id: Series[pd.StringDtype] = pa.Field(unique=True)
    model_id: Series[pd.StringDtype]
    cell_line_name: Series[pd.StringDtype]
    stripped_cell_line_name: Series[pd.StringDtype]
    model_condition_id: Series[pd.StringDtype]
    omics_profile_id: Series[pd.StringDtype]
    delivery_cram_bam: Series[pd.StringDtype] = pa.Field(nullable=True)
    delivery_crai_bai: Series[pd.StringDtype] = pa.Field(nullable=True)
    delivery_file_format: Series[pd.StringDtype] = pa.Field(
        isin={"CRAM", "BAM"}, nullable=True
    )
    delivery_cram_bam_size: Series[pd.Int64Dtype] = pa.Field(unique=True)
    delivery_ref: Series[pd.StringDtype]
    delivery_ref_alt: Series[pd.StringDtype]
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
    analysis_ready_bam_size: Series[pd.Int64Dtype] = pa.Field(nullable=True)
    ref: Series[pd.StringDtype]
    ref_alt: Series[pd.StringDtype]
    ref_amb: Series[pd.StringDtype]
    ref_ann: Series[pd.StringDtype]
    ref_bwt: Series[pd.StringDtype]
    ref_dict: Series[pd.StringDtype]
    ref_fasta: Series[pd.StringDtype]
    ref_fasta_index: Series[pd.StringDtype]
    ref_pac: Series[pd.StringDtype]
    ref_sa: Series[pd.StringDtype]
    automation_status: Series[pd.StringDtype] = pa.Field(nullable=True)


class GcsObject(CoercedDataFrame):
    url: Series[pd.StringDtype]
    size: Series[pd.Int64Dtype]
    crc32c: Series[pd.StringDtype]
    gcs_obj_updated_at: Series[pd.Timestamp]


class CopiedSampleFiles(CoercedDataFrame):
    sample_id: Series[pd.StringDtype]
    url_kind: Series[pd.StringDtype]
    new_url: Series[pd.StringDtype]
    url: Series[pd.StringDtype]
    copied: Series[pd.BooleanDtype]


class AlignedSamples(CoercedDataFrame):
    sample_id: Series[pd.StringDtype] = pa.Field(unique=True)
    analysis_ready_bam: Series[pd.StringDtype] = pa.Field(unique=True)
    analysis_ready_bai: Series[pd.StringDtype] = pa.Field(unique=True)


class AlignedSamplesWithObjectMetadata(AlignedSamples):
    crc32c: Series[pd.StringDtype] = pa.Field(unique=True)
    size: Series[pd.Int64Dtype] = pa.Field(unique=True)


class NewSequencingAlignments(CoercedDataFrame):
    omics_sequencing_id: Series[pd.StringDtype] = pa.Field(unique=True)
    url: Series[pd.StringDtype] = pa.Field(unique=True)
    index_url: Series[pd.StringDtype] = pa.Field(unique=True)
    sequencing_alignment_source: Series[pd.StringDtype]
    reference_genome: Series[pd.StringDtype]
    crc32c_hash: Series[pd.StringDtype] = pa.Field(unique=True)
    size: Series[pd.Int64Dtype] = pa.Field(unique=True)


class DeltaJob(BaseModel):
    workflow_name: str
    entity_type: str
    entity_set_type: str
    entity_id_col: str
    expression: str
    input_cols: set[str] | None = None
    output_cols: set[str] | None = None
    resubmit_n_times: int = 0
    force_retry: bool = False
    use_callcache: bool = True
    use_reference_disks: bool = False
    memory_retry_multiplier: float = 1.0
    max_n_entities: int | None = None
    dry_run: bool = False


PydanticBaseModel = TypeVar("PydanticBaseModel", bound=BaseModel)
PanderaBaseSchema = TypeVar("PanderaBaseSchema", bound=CoercedDataFrame)
TypedDataFrame = pandera.typing.DataFrame
T = TypeVar("T")
