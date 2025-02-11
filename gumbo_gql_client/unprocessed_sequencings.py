from typing import List

from .base_model import BaseModel


class UnprocessedSequencings(BaseModel):
    records: List["UnprocessedSequencingsRecords"]


class UnprocessedSequencingsRecords(BaseModel):
    omics_sequencing_id: str


UnprocessedSequencings.model_rebuild()
