import json
import tempfile
from pathlib import Path

from click import echo
from firecloud import api as firecloud_api
from google.cloud import storage

from omics_wgs_pipeline.terra import call_firecloud_api
from omics_wgs_pipeline.wdl import persist_wdl_script


class TerraWorkflow:
    def __init__(
        self,
        gcp_project_id: str,
        pipelines_bucket_name: str,
        repo_namespace: str,
        repo_method_name: str,
        workspace_namespace: str,
        workspace_name: str,
        method_config_name: str,
        method_synopsis: str,
        workflow_wdl_path: Path,
        method_config_json_path: Path,
        firecloud_owners: list[str],
    ) -> None:
        self.gcp_project_id = gcp_project_id
        self.pipelines_bucket_name = pipelines_bucket_name
        self.repo_namespace = repo_namespace
        self.repo_method_name = repo_method_name
        self.workspace_namespace = workspace_namespace
        self.workspace_name = workspace_name
        self.method_config_name = method_config_name
        self.method_synopsis = method_synopsis
        self.workflow_wdl_path = workflow_wdl_path
        self.method_config_json_path = method_config_json_path
        self.firecloud_owners = firecloud_owners
        self.method_config = json.load(open(self.method_config_json_path, "r"))

    def persist_workflow(self) -> str:
        storage_client = storage.Client(project=self.gcp_project_id)
        bucket = storage_client.bucket(
            self.pipelines_bucket_name, user_project=self.gcp_project_id
        )
        wdl = persist_wdl_script(
            bucket=bucket,
            wdl_path=self.workflow_wdl_path,
            subpath="wdl",
        )["wdl"]

        return wdl

    def create_method(self) -> dict:
        """
        Create the initial method using the WDL file in this repo.

        :return: the latest method's snapshot
        """

        snapshot = call_firecloud_api(
            firecloud_api.update_repository_method,
            namespace=self.repo_namespace,
            method=self.repo_method_name,
            synopsis=self.method_synopsis,
            wdl=self.workflow_wdl_path,
        )

        # set permissions
        call_firecloud_api(
            firecloud_api.update_repository_method_acl,
            namespace=self.repo_namespace,
            method=self.repo_method_name,
            snapshot_id=snapshot["snapshotId"],
            acl_updates=[{"user": x, "role": "OWNER"} for x in self.firecloud_owners],
        )

        return snapshot

    def update_method(self) -> dict:
        """
        Update the Terra method using the WDL file in this repo.

        :return: the latest method's snapshot
        """

        wdl = self.persist_workflow()

        with tempfile.NamedTemporaryFile("w") as f:
            f.write(wdl)

            snapshot = call_firecloud_api(
                firecloud_api.update_repository_method,
                namespace=self.repo_namespace,
                method=self.repo_method_name,
                synopsis=self.method_synopsis,
                wdl=f.name,
            )

        # set permissions again
        call_firecloud_api(
            firecloud_api.update_repository_method_acl,
            namespace=self.repo_namespace,
            method=self.repo_method_name,
            snapshot_id=snapshot["snapshotId"],
            acl_updates=[{"user": x, "role": "OWNER"} for x in self.firecloud_owners],
        )

        return snapshot

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

    def update_workspace_config(self, config_body: dict) -> None:
        """
        Update the Terra workspace config for a given method.

        :param config_body: a dictionary containing the method config
        """

        call_firecloud_api(
            firecloud_api.update_workspace_config,
            namespace=self.workspace_namespace,
            workspace=self.workspace_name,
            cnamespace=self.repo_namespace,
            configname=self.repo_method_name,
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

    def update_workflow(self) -> None:
        """
        Update the Terra workflow (method and method config).
        """

        snapshots = self.get_method_snapshots()

        # update or create the method for the current WDL file
        if len(snapshots) == 0:
            snapshot = self.create_method()
        else:
            snapshot = self.update_method()

        # assocate the method config with the latest method version
        self.method_config["methodRepoMethod"]["methodVersion"] = snapshot["snapshotId"]

        # check if there is already a workspace config
        res = firecloud_api.get_workspace_config(
            namespace=self.workspace_namespace,
            workspace=self.workspace_name,
            cnamespace=self.repo_namespace,
            config=self.method_config_name,
        )

        # update or create the method config
        if res.status_code == 404:
            self.create_workspace_config(self.method_config)
        else:
            self.update_workspace_config(self.method_config)

        # don't let old method configs accumulate
        self.delete_old_method_snapshots(2)
