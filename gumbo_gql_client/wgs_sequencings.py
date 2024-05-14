from typing import List, Optional

from pydantic import Field

from .base_model import BaseModel


class WgsSequencings(BaseModel):
    records: List["WgsSequencingsRecords"]


class WgsSequencingsRecords(BaseModel):
    hg_19_bai_filepath: Optional[str] = Field(alias="hg19_bai_filepath")
    hg_19_bam_filepath: Optional[str] = Field(alias="hg19_bam_filepath")
    hg_38_crai_filepath: Optional[str] = Field(alias="hg38_crai_filepath")
    hg_38_cram_filepath: Optional[str] = Field(alias="hg38_cram_filepath")
    bai_filepath: Optional[str]
    bam_filepath: Optional[str]
    sequencing_id: str


WgsSequencings.model_rebuild()
