from typing import List, Optional

from .base_model import BaseModel


class UpdateSequencingAlignment(BaseModel):
    set_username: List["UpdateSequencingAlignmentSetUsername"]
    update_sequencing_alignment_by_pk: Optional[
        "UpdateSequencingAlignmentUpdateSequencingAlignmentByPk"
    ]


class UpdateSequencingAlignmentSetUsername(BaseModel):
    username: str


class UpdateSequencingAlignmentUpdateSequencingAlignmentByPk(BaseModel):
    id: int


UpdateSequencingAlignment.model_rebuild()
