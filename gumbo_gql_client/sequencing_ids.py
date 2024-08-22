from typing import List

from .base_model import BaseModel


class SequencingIds(BaseModel):
    records: List["SequencingIdsRecords"]


class SequencingIdsRecords(BaseModel):
    sequencing_id: str


SequencingIds.model_rebuild()
