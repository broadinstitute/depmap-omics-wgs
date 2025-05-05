from typing import List, Optional

from .base_model import BaseModel


class WgsSequencingAlignments(BaseModel):
    records: List["WgsSequencingAlignmentsRecords"]


class WgsSequencingAlignmentsRecords(BaseModel):
    index_url: Optional[str]
    reference_genome: Optional[str]
    omics_sequencing_id: str
    sequencing_alignment_source: str
    url: str
    size: int


WgsSequencingAlignments.model_rebuild()
