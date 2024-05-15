from typing import List, Optional

from .base_model import BaseModel


class SequencingTaskEntities(BaseModel):
    records: List["SequencingTaskEntitiesRecords"]


class SequencingTaskEntitiesRecords(BaseModel):
    id: int
    sequencing_id: Optional[str]


SequencingTaskEntities.model_rebuild()
