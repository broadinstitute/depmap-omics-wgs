from typing import TypedDict, TypeVar

import pandera.typing
from pydantic import BaseModel

from omics_wgs_pipeline.validators import CoercedDataFrame

PydanticBaseModel = TypeVar("PydanticBaseModel", bound=BaseModel)
PanderaBaseSchema = TypeVar("PanderaBaseSchema", bound=CoercedDataFrame)
TypedDataFrame = pandera.typing.DataFrame


class PersistedWdl(TypedDict):
    wdl: str
    public_url: str


class TerraJobSubmissionKwargs(TypedDict):
    entity: str | None
    etype: str | None
    expression: str | None
    use_callcache: bool | None
    delete_intermediate_output_files: bool | None
    use_reference_disks: bool | None
    memory_retry_multiplier: float | None
    workflow_failure_mode: str | None
    user_comment: str | None
