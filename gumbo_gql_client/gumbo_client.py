from typing import Any, Dict, List

from .base_client import BaseClient
from .get_task_results import GetTaskResults
from .input_types import (
    sequencing_alignment_insert_input,
    task_entity_insert_input,
    task_result_bool_exp,
    task_result_insert_input,
)
from .insert_sequencing_alignments import InsertSequencingAlignments
from .insert_task_entities import InsertTaskEntities
from .insert_task_results import InsertTaskResults
from .insert_terra_sync import InsertTerraSync
from .sequencing_ids import SequencingIds
from .sequencing_task_entities import SequencingTaskEntities
from .unprocessed_sequencings import UnprocessedSequencings
from .wgs_sequencing_alignments import WgsSequencingAlignments


def gql(q: str) -> str:
    return q


class GumboClient(BaseClient):
    def sequencing_ids(self, **kwargs: Any) -> SequencingIds:
        query = gql(
            """
            query SequencingIds {
              records: omics_sequencing {
                sequencing_id
              }
            }
            """
        )
        variables: Dict[str, object] = {}
        response = self.execute(
            query=query, operation_name="SequencingIds", variables=variables, **kwargs
        )
        data = self.get_data(response)
        return SequencingIds.model_validate(data)

    def unprocessed_sequencings(self, **kwargs: Any) -> UnprocessedSequencings:
        query = gql(
            """
            query UnprocessedSequencings {
              records: omics_sequencing(
                where: {_and: [{expected_type: {_eq: "wgs"}, blacklist: {_neq: true}, _or: [{processed_sequence: {_neq: true}}, {processed_sequence: {_is_null: true}}], omics_profile: {blacklist_omics: {_neq: true}}}]}
              ) {
                omics_sequencing_id: sequencing_id
              }
            }
            """
        )
        variables: Dict[str, object] = {}
        response = self.execute(
            query=query,
            operation_name="UnprocessedSequencings",
            variables=variables,
            **kwargs
        )
        data = self.get_data(response)
        return UnprocessedSequencings.model_validate(data)

    def wgs_sequencing_alignments(self, **kwargs: Any) -> WgsSequencingAlignments:
        query = gql(
            """
            query WgsSequencingAlignments {
              records: sequencing_alignment(
                where: {_and: [{omics_sequencing: {expected_type: {_eq: "wgs"}, blacklist: {_neq: true}, omics_profile: {blacklist_omics: {_neq: true}}}}]}
              ) {
                index_url
                reference_genome
                omics_sequencing_id
                sequencing_alignment_source
                url
                size
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

    def sequencing_task_entities(self, **kwargs: Any) -> SequencingTaskEntities:
        query = gql(
            """
            query SequencingTaskEntities {
              records: task_entity(where: {sequencing_id: {_is_null: false}}) {
                id
                sequencing_id
              }
            }
            """
        )
        variables: Dict[str, object] = {}
        response = self.execute(
            query=query,
            operation_name="SequencingTaskEntities",
            variables=variables,
            **kwargs
        )
        data = self.get_data(response)
        return SequencingTaskEntities.model_validate(data)

    def insert_task_entities(
        self, username: str, objects: List[task_entity_insert_input], **kwargs: Any
    ) -> InsertTaskEntities:
        query = gql(
            """
            mutation InsertTaskEntities($username: String!, $objects: [task_entity_insert_input!]!) {
              set_username(args: {_username: $username}) {
                username
              }
              insert_task_entity(objects: $objects) {
                returning {
                  id
                  sequencing_id
                }
              }
            }
            """
        )
        variables: Dict[str, object] = {"username": username, "objects": objects}
        response = self.execute(
            query=query,
            operation_name="InsertTaskEntities",
            variables=variables,
            **kwargs
        )
        data = self.get_data(response)
        return InsertTaskEntities.model_validate(data)

    def get_task_results(
        self, where: task_result_bool_exp, **kwargs: Any
    ) -> GetTaskResults:
        query = gql(
            """
            query GetTaskResults($where: task_result_bool_exp!) {
              records: task_result(where: $where) {
                id
                crc32c_hash
                completed_at
                format
                label
                size
                task_entity {
                  id
                  sequencing_id
                }
                terra_entity_name
                terra_entity_type
                terra_method_config_name
                terra_method_config_namespace
                terra_submission_id
                terra_workflow_id
                terra_workflow_inputs
                terra_workflow_root_dir
                terra_workspace_id
                terra_workspace_name
                terra_workspace_namespace
                url
                value
                workflow_name
                workflow_source_url
                workflow_version
              }
            }
            """
        )
        variables: Dict[str, object] = {"where": where}
        response = self.execute(
            query=query, operation_name="GetTaskResults", variables=variables, **kwargs
        )
        data = self.get_data(response)
        return GetTaskResults.model_validate(data)

    def insert_task_results(
        self, username: str, objects: List[task_result_insert_input], **kwargs: Any
    ) -> InsertTaskResults:
        query = gql(
            """
            mutation InsertTaskResults($username: String!, $objects: [task_result_insert_input!]!) {
              set_username(args: {_username: $username}) {
                username
              }
              insert_task_result(
                objects: $objects
                on_conflict: {constraint: task_result_pkey, update_columns: []}
              ) {
                affected_rows
              }
            }
            """
        )
        variables: Dict[str, object] = {"username": username, "objects": objects}
        response = self.execute(
            query=query,
            operation_name="InsertTaskResults",
            variables=variables,
            **kwargs
        )
        data = self.get_data(response)
        return InsertTaskResults.model_validate(data)

    def insert_terra_sync(
        self,
        username: str,
        terra_workspace_namespace: str,
        terra_workspace_name: str,
        task_results: List[task_result_insert_input],
        **kwargs: Any
    ) -> InsertTerraSync:
        query = gql(
            """
            mutation InsertTerraSync($username: String!, $terra_workspace_namespace: String!, $terra_workspace_name: String!, $task_results: [task_result_insert_input!]!) {
              set_username(args: {_username: $username}) {
                username
              }
              insert_terra_sync(
                objects: {terra_workspace_name: $terra_workspace_name, terra_workspace_namespace: $terra_workspace_namespace, task_results: {data: $task_results, on_conflict: {constraint: task_result_pkey, update_columns: []}}}
              ) {
                returning {
                  id
                }
              }
            }
            """
        )
        variables: Dict[str, object] = {
            "username": username,
            "terra_workspace_namespace": terra_workspace_namespace,
            "terra_workspace_name": terra_workspace_name,
            "task_results": task_results,
        }
        response = self.execute(
            query=query, operation_name="InsertTerraSync", variables=variables, **kwargs
        )
        data = self.get_data(response)
        return InsertTerraSync.model_validate(data)
