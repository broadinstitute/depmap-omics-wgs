import datetime
import json
import tempfile
from pathlib import Path
from typing import Iterable, Type, TypeVar, Unpack

import pandas as pd
from click import echo
from firecloud import api as firecloud_api
from google.cloud import storage

from omics_wgs_pipeline.types import TerraJobSubmissionKwargs, TypedDataFrame
from omics_wgs_pipeline.utils import call_firecloud_api
from omics_wgs_pipeline.validators import CoercedDataFrame
from omics_wgs_pipeline.wdl import persist_wdl_script

T = TypeVar("T", bound=CoercedDataFrame)


class TerraWorkflow:
    def __init__(
        self,
        gcp_project_id: str,
        pipelines_bucket_name: str,
        repo_namespace: str,
        repo_method_name: str,
        method_config_name: str,
        method_synopsis: str,
        workflow_wdl_path: Path,
        method_config_json_path: Path,
    ) -> None:
        self.gcp_project_id = gcp_project_id
        self.pipelines_bucket_name = pipelines_bucket_name
        self.repo_namespace = repo_namespace
        self.repo_method_name = repo_method_name
        self.method_config_name = method_config_name
        self.method_synopsis = method_synopsis
        self.workflow_wdl_path = workflow_wdl_path
        self.method_config_json_path = method_config_json_path
        self.method_config = json.load(open(self.method_config_json_path, "r"))
        self.persisted_wdl_url: str | None = None
        self.persisted_wdl: str | None = None

    def persist_method_on_gcs(self) -> None:
        """
        Upload the method's WDL script to GCS, rewriting import statements for dependent
        WDL scripts as needed.
        """

        storage_client = storage.Client(project=self.gcp_project_id)
        bucket = storage_client.bucket(
            self.pipelines_bucket_name, user_project=self.gcp_project_id
        )

        res = persist_wdl_script(
            bucket=bucket, wdl_path=self.workflow_wdl_path, subpath="wdl"
        )

        self.persisted_wdl = res["wdl"]
        self.persisted_wdl_url = res["public_url"]

    def get_method_snapshots(self) -> list[dict]:
        """
        Get all of the snapshots of the method.

        :return: list of snapshot information, most recent first
        """

        snapshots = call_firecloud_api(
            firecloud_api.list_repository_methods,
            namespace=self.repo_namespace,
            name=self.repo_method_name,
        )

        snapshots.sort(key=lambda x: x["snapshotId"], reverse=True)
        return snapshots

    def delete_old_method_snapshots(self, n: int = 2) -> None:
        """
        Delete all but `n` of the most recent snapshots of the method. This might fail
        if the service account doesn't have OWNER permission on the method configuration
        namespace. This can't be set using the firecloud package, so POST to the
        `setConfigNamespaceACL` endpoint on https://api.firecloud.org/ once to resolve
        this problem.

        :param n: the number of snapshots to keep
        """

        snapshots = self.get_method_snapshots()

        to_delete = snapshots[n:]
        echo(f"Deleting {len(to_delete)} old snapshot(s)")

        for x in to_delete:
            call_firecloud_api(
                firecloud_api.delete_repository_method,
                namespace=self.repo_namespace,
                name=self.repo_method_name,
                snapshot_id=x["snapshotId"],
            )


class TerraWorkspace:
    def __init__(
        self, workspace_namespace: str, workspace_name: str, firecloud_owners: list[str]
    ) -> None:
        self.workspace_namespace = workspace_namespace
        self.workspace_name = workspace_name
        self.firecloud_owners = firecloud_owners

    def get_entities(
        self, entity_type: str, pandera_schema: Type[T]
    ) -> TypedDataFrame[T]:
        """
        Get a data frame of entities from a Terra data table.

        :param entity_type: the kind of entity (e.g. "sample")
        :param pandera_schema: a Pandera schema for the output data frame
        """

        j = call_firecloud_api(
            firecloud_api.get_entities,
            namespace=self.workspace_namespace,
            workspace=self.workspace_name,
            etype=entity_type,
        )

        return TypedDataFrame[pandera_schema](pd.DataFrame(j))

    def upload_entities(self, df: pd.DataFrame) -> None:
        """
        Upload a data frame of entities to a Terra data table.

        :param df: a data frame of entities
        """

        with tempfile.NamedTemporaryFile(suffix="tsv") as f:
            df.to_csv(f, sep="\t", index=False)

            echo(f"Upserting {len(df)} entities to Terra")

            call_firecloud_api(
                firecloud_api.upload_entities_tsv,
                namespace=self.workspace_namespace,
                workspace=self.workspace_name,
                entities_tsv=f.name,
                model="flexible",
            )

    def create_workspace_config(self, config_body: dict) -> None:
        """
        Create a Terra workspace config for the pipeline in Terra.

        :param config_body: a dictionary containing the method config
        """

        call_firecloud_api(
            firecloud_api.create_workspace_config,
            namespace=self.workspace_namespace,
            workspace=self.workspace_name,
            body=config_body,
        )

        # set permissions
        call_firecloud_api(
            firecloud_api.update_workspace_acl,
            namespace=self.workspace_namespace,
            workspace=self.workspace_name,
            acl_updates=[
                {"email": x, "accessLevel": "OWNER"} for x in self.firecloud_owners
            ],
        )

    def update_workspace_config(
        self, terra_workflow: TerraWorkflow, config_body: dict
    ) -> None:
        """
        Update the Terra workspace config for a given method.

        :param terra_workflow: a `TerraWorkflow` instance
        :param config_body: a dictionary containing the method config
        """

        call_firecloud_api(
            firecloud_api.update_workspace_config,
            namespace=self.workspace_namespace,
            workspace=self.workspace_name,
            cnamespace=terra_workflow.repo_namespace,
            configname=terra_workflow.repo_method_name,
            body=config_body,
        )

        # set permissions again
        call_firecloud_api(
            firecloud_api.update_workspace_acl,
            namespace=self.workspace_namespace,
            workspace=self.workspace_name,
            acl_updates=[
                {"email": x, "accessLevel": "OWNER"} for x in self.firecloud_owners
            ],
        )

    def create_method(self, terra_workflow: TerraWorkflow) -> dict:
        """
        Create the initial method using the WDL file in this repo.

        :param terra_workflow: a `TerraWorkflow` instance
        :return: the latest method's snapshot
        """

        snapshot = call_firecloud_api(
            firecloud_api.update_repository_method,
            namespace=terra_workflow.repo_namespace,
            method=terra_workflow.repo_method_name,
            synopsis=terra_workflow.method_synopsis,
            wdl=terra_workflow.workflow_wdl_path,
        )

        # set permissions
        call_firecloud_api(
            firecloud_api.update_repository_method_acl,
            namespace=terra_workflow.repo_namespace,
            method=terra_workflow.repo_method_name,
            snapshot_id=snapshot["snapshotId"],
            acl_updates=[{"user": x, "role": "OWNER"} for x in self.firecloud_owners],
        )

        return snapshot

    def update_method(self, terra_workflow: TerraWorkflow) -> dict:
        """
        Update the Terra method using a WDL file in this repo.

        :param terra_workflow: a `TerraWorkflow` instance
        :return: the latest method's snapshot
        """

        # get contents of WDL uploaded to GCS
        terra_workflow.persist_method_on_gcs()
        assert terra_workflow.persisted_wdl is not None

        with tempfile.NamedTemporaryFile("w") as f:
            f.write(terra_workflow.persisted_wdl)

            snapshot = call_firecloud_api(
                firecloud_api.update_repository_method,
                namespace=terra_workflow.repo_namespace,
                method=terra_workflow.repo_method_name,
                synopsis=terra_workflow.method_synopsis,
                wdl=f.name,
            )

        # set permissions again
        call_firecloud_api(
            firecloud_api.update_repository_method_acl,
            namespace=terra_workflow.repo_namespace,
            method=terra_workflow.repo_method_name,
            snapshot_id=snapshot["snapshotId"],
            acl_updates=[{"user": x, "role": "OWNER"} for x in self.firecloud_owners],
        )

        return snapshot

    def update_workflow(self, terra_workflow: TerraWorkflow) -> None:
        """
        Update the Terra workflow (method and method config).

        :param terra_workflow: a `TerraWorkflow` instance
        """

        snapshots = terra_workflow.get_method_snapshots()

        # update or create the method for the current WDL file
        if len(snapshots) == 0:
            snapshot = self.create_method(terra_workflow)
        else:
            snapshot = self.update_method(terra_workflow)

        # assocate the method config with the latest method version
        terra_workflow.method_config["methodRepoMethod"]["methodVersion"] = snapshot[
            "snapshotId"
        ]

        # check if there is already a workspace config
        res = firecloud_api.get_workspace_config(
            namespace=self.workspace_namespace,
            workspace=self.workspace_name,
            cnamespace=terra_workflow.repo_namespace,
            config=terra_workflow.method_config_name,
        )

        # update or create the method config
        if res.status_code == 404:
            self.create_workspace_config(terra_workflow.method_config)
        else:
            self.update_workspace_config(terra_workflow, terra_workflow.method_config)

        # don't let old method configs accumulate
        terra_workflow.delete_old_method_snapshots(2)

    def submit_workflow_run(
        self, terra_workflow: TerraWorkflow, **kwargs: Unpack[TerraJobSubmissionKwargs]
    ) -> None:
        """
        Submit a run of a workflow.

        :param terra_workflow: a `TerraWorkflow` instance
        """

        call_firecloud_api(
            firecloud_api.create_submission,
            wnamespace=self.workspace_namespace,
            workspace=self.workspace_name,
            cnamespace=terra_workflow.repo_namespace,
            config=terra_workflow.repo_method_name,
            **kwargs,
        )

    def create_sample_set(self, sample_ids: Iterable[str], suffix: str) -> str:
        """
        Create a new sample set for a list of sample IDs and upload it to Terra.

        :param sample_ids: a list of sample IDs
        :param suffix: a suffix to add to the sample set ID (e.g.
        "preprocess_wgs_sample")
        :return: the ID of the new sample set
        """

        # make an ID for the sample set of new samples
        sample_set_id = "_".join(
            [
                "samples",
                datetime.datetime.now(datetime.UTC)
                .isoformat(timespec="seconds")
                .replace(":", "-"),
                suffix,
            ]
        )

        # construct a data frame of sample IDs for this sample set
        sample_sets = pd.DataFrame({"entity:sample_id": sample_ids}, dtype="string")
        sample_sets["entity:sample_set_id"] = sample_set_id

        echo("Creating new sample set in Terra")
        self.upload_entities(
            sample_sets.loc[:, ["entity:sample_set_id"]].drop_duplicates()
        )

        # construct the join table between the sample set and its samples
        sample_sets = sample_sets.rename(
            columns={
                "entity:sample_set_id": "membership:sample_set_id",
                "entity:sample_id": "sample",
            }
        )

        sample_sets = sample_sets.loc[:, ["membership:sample_set_id", "sample"]]

        echo(f"Adding {len(sample_sets)} samples to sample set {sample_set_id}")
        self.upload_entities(sample_sets)

        return sample_set_id
