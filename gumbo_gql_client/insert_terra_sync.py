from typing import List, Optional

from .base_model import BaseModel


class InsertTerraSync(BaseModel):
    set_username: List["InsertTerraSyncSetUsername"]
    insert_terra_sync: Optional["InsertTerraSyncInsertTerraSync"]


class InsertTerraSyncSetUsername(BaseModel):
    username: str


class InsertTerraSyncInsertTerraSync(BaseModel):
    returning: List["InsertTerraSyncInsertTerraSyncReturning"]


class InsertTerraSyncInsertTerraSyncReturning(BaseModel):
    id: int


InsertTerraSync.model_rebuild()
InsertTerraSyncInsertTerraSync.model_rebuild()
