from typing import List, Optional

from .base_model import BaseModel


class WgsSequencingAlignments(BaseModel):
    records: List["WgsSequencingAlignmentsRecords"]


class WgsSequencingAlignmentsRecords(BaseModel):
    model_id: Optional[str]
    model_condition_id: Optional[str]
    omics_profile_id: Optional[str]
    omics_sequencing_id: Optional[str]
    model: Optional["WgsSequencingAlignmentsRecordsModel"]
    omics_sequencing: Optional["WgsSequencingAlignmentsRecordsOmicsSequencing"]


class WgsSequencingAlignmentsRecordsModel(BaseModel):
    cell_line_name: Optional[str]
    stripped_cell_line_name: Optional[str]


class WgsSequencingAlignmentsRecordsOmicsSequencing(BaseModel):
    sequencing_alignments: List[
        "WgsSequencingAlignmentsRecordsOmicsSequencingSequencingAlignments"
    ]


class WgsSequencingAlignmentsRecordsOmicsSequencingSequencingAlignments(BaseModel):
    sequencing_alignment_id: int
    url: str
    index_url: Optional[str]
    reference_genome: Optional[str]
    sequencing_alignment_source: str
    size: int


WgsSequencingAlignments.model_rebuild()
WgsSequencingAlignmentsRecords.model_rebuild()
WgsSequencingAlignmentsRecordsOmicsSequencing.model_rebuild()
