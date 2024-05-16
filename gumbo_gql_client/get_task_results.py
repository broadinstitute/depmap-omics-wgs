from typing import Any, List, Optional

from pydantic import Field

from .base_model import BaseModel


class GetTaskResults(BaseModel):
    records: List["GetTaskResultsRecords"]


class GetTaskResultsRecords(BaseModel):
    id: int
    crc_32_c_hash: Optional[str] = Field(alias="crc32c_hash")
    created_at: Any
    format: Optional[str]
    label: str
    size: Optional[int]
    task_entity: "GetTaskResultsRecordsTaskEntity"
    task_name: str
    terra_entity_name: Optional[str]
    terra_entity_type: Optional[str]
    terra_method_config_name: Optional[str]
    terra_method_config_namespace: Optional[str]
    terra_submission_id: Optional[str]
    terra_workflow_id: Optional[str]
    terra_workflow_inputs: Optional[Any]
    terra_workflow_root_dir: Optional[str]
    terra_workspace_id: Optional[str]
    terra_workspace_name: Optional[str]
    terra_workspace_namespace: Optional[str]
    url: Optional[str]
    workflow_name: str
    workflow_source_url: Optional[str]


class GetTaskResultsRecordsTaskEntity(BaseModel):
    id: int
    sequencing_id: Optional[str]


GetTaskResults.model_rebuild()
GetTaskResultsRecords.model_rebuild()
