from typing import Any, Dict

from .base_client import BaseClient
from .wgs_sequencings import WgsSequencings


def gql(q: str) -> str:
    return q


class GumboClient(BaseClient):
    def wgs_sequencings(self, **kwargs: Any) -> WgsSequencings:
        query = gql(
            """
            query WgsSequencings {
              records: omics_sequencing(
                where: {_and: [{_or: [{blacklist: {_is_null: true}}, {blacklist: {_neq: true}}]}, {_or: [{bam_filepath: {_is_null: false}}, {hg19_bam_filepath: {_is_null: false}}, {hg38_cram_filepath: {_is_null: false}}]}, {expected_type: {_eq: "wgs"}}, {omics_profile: {_or: [{blacklist_omics: {_is_null: true}}, {blacklist_omics: {_neq: true}}]}}]}
              ) {
                hg19_bai_filepath
                hg19_bam_filepath
                hg38_crai_filepath
                hg38_cram_filepath
                bai_filepath
                bam_filepath
                sequencing_id
              }
            }
            """
        )
        variables: Dict[str, object] = {}
        response = self.execute(
            query=query, operation_name="WgsSequencings", variables=variables, **kwargs
        )
        data = self.get_data(response)
        return WgsSequencings.model_validate(data)
