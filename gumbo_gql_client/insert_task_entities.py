from typing import List, Optional

from .base_model import BaseModel


class InsertTaskEntities(BaseModel):
    set_username: List["InsertTaskEntitiesSetUsername"]
    insert_task_entity: Optional["InsertTaskEntitiesInsertTaskEntity"]


class InsertTaskEntitiesSetUsername(BaseModel):
    username: str


class InsertTaskEntitiesInsertTaskEntity(BaseModel):
    returning: List["InsertTaskEntitiesInsertTaskEntityReturning"]


class InsertTaskEntitiesInsertTaskEntityReturning(BaseModel):
    id: int
    omics_sequencing_id: Optional[str]


InsertTaskEntities.model_rebuild()
InsertTaskEntitiesInsertTaskEntity.model_rebuild()
