from typing import List, Optional

from .base_model import BaseModel


class InsertTaskResults(BaseModel):
    set_username: List["InsertTaskResultsSetUsername"]
    insert_task_result: Optional["InsertTaskResultsInsertTaskResult"]


class InsertTaskResultsSetUsername(BaseModel):
    username: str


class InsertTaskResultsInsertTaskResult(BaseModel):
    affected_rows: int


InsertTaskResults.model_rebuild()
