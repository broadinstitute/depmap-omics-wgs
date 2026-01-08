from typing import Any, Dict, List

from .base_client import BaseClient
from .input_types import (
    sequencing_alignment_insert_input,
    sequencing_alignment_set_input,
)
from .insert_sequencing_alignments import InsertSequencingAlignments
from .update_sequencing_alignment import UpdateSequencingAlignment
from .wgs_sequencing_alignments import WgsSequencingAlignments


def gql(q: str) -> str:
    return q


class GumboClient(BaseClient):
    def wgs_sequencing_alignments(self, **kwargs: Any) -> WgsSequencingAlignments:
        query = gql(
            """
            query WgsSequencingAlignments {
              records: omics_mapping(where: {datatype: {_eq: "wgs"}}) {
                model_id
                model_condition_id
                omics_profile_id
                omics_sequencing_id
                model {
                  cell_line_name
                  stripped_cell_line_name
                }
                omics_sequencing {
                  sequencing_alignments {
                    sequencing_alignment_id: id
                    url
                    index_url
                    reference_genome
                    sequencing_alignment_source
                    size
                  }
                }
              }
            }
            """
        )
        variables: Dict[str, object] = {}
        response = self.execute(
            query=query,
            operation_name="WgsSequencingAlignments",
            variables=variables,
            **kwargs
        )
        data = self.get_data(response)
        return WgsSequencingAlignments.model_validate(data)

    def insert_sequencing_alignments(
        self,
        username: str,
        objects: List[sequencing_alignment_insert_input],
        **kwargs: Any
    ) -> InsertSequencingAlignments:
        query = gql(
            """
            mutation InsertSequencingAlignments($_username: String!, $objects: [sequencing_alignment_insert_input!]!) {
              set_username(args: {_username: $_username}) {
                username
              }
              insert_sequencing_alignment(objects: $objects) {
                affected_rows
              }
            }
            """
        )
        variables: Dict[str, object] = {"_username": username, "objects": objects}
        response = self.execute(
            query=query,
            operation_name="InsertSequencingAlignments",
            variables=variables,
            **kwargs
        )
        data = self.get_data(response)
        return InsertSequencingAlignments.model_validate(data)

    def update_sequencing_alignment(
        self,
        username: str,
        id: int,
        object: sequencing_alignment_set_input,
        **kwargs: Any
    ) -> UpdateSequencingAlignment:
        query = gql(
            """
            mutation UpdateSequencingAlignment($_username: String!, $id: bigint!, $object: sequencing_alignment_set_input!) {
              set_username(args: {_username: $_username}) {
                username
              }
              update_sequencing_alignment_by_pk(pk_columns: {id: $id}, _set: $object) {
                id
              }
            }
            """
        )
        variables: Dict[str, object] = {
            "_username": username,
            "id": id,
            "object": object,
        }
        response = self.execute(
            query=query,
            operation_name="UpdateSequencingAlignment",
            variables=variables,
            **kwargs
        )
        data = self.get_data(response)
        return UpdateSequencingAlignment.model_validate(data)
