import datetime
import json
import pathlib
import tempfile
from pathlib import Path
from typing import Any, Callable, Iterable, Type, Unpack

import pandas as pd
import requests
from click import echo
from firecloud import api as firecloud_api

from gumbo_gql_client import task_result_insert_input
from omics_wgs_pipeline.types import (
    PanderaBaseSchema,
    PersistedWdl,
    TerraJobSubmissionKwargs,
    TypedDataFrame,
)
from omics_wgs_pipeline.utils import batch_evenly, expand_dict_columns, maybe_retry
from omics_wgs_pipeline.wdl import GistedWdl


class TerraWorkflow:
    def __init__(
        self,
        github_pat: str,
        pipelines_bucket_name: str,
        repo_namespace: str,
        repo_method_name: str,
        method_config_name: str,
        method_synopsis: str,
        workflow_wdl_path: Path,
        method_config_json_path: Path,
    ) -> None:
        self.github_pat = github_pat
        self.pipelines_bucket_name = pipelines_bucket_name
        self.repo_namespace = repo_namespace
        self.repo_method_name = repo_method_name
        self.method_config_name = method_config_name
        self.method_synopsis = method_synopsis
        self.workflow_wdl_path = workflow_wdl_path
        self.method_config_json_path = method_config_json_path
        self.method_config = json.load(open(self.method_config_json_path, "r"))
        self.persisted_wdl_script: PersistedWdl | None = None

    def persist_method_on_gcs(self) -> None:
        """
        Upload the method's WDL script to GitHub, rewriting import statements for
        dependent WDL scripts as needed.
        """

        if self.persisted_wdl_script is None:
            echo(f"Persisting {self.workflow_wdl_path} on GitHub")
            gisted_wdl = GistedWdl(
                method_name=self.repo_method_name, github_pat=self.github_pat
            )
            self.persisted_wdl_script = gisted_wdl.persist_wdl_script(
                wdl_path=self.workflow_wdl_path, subpath="wdl"
            )

    def get_method_snapshots(self) -> list[dict]:
        """
        Get all of the snapshots of the method.

        :return: list of snapshot information, most recent first
        """

        echo(f"Getting {self.repo_method_name} method snapshots")
        snapshots = call_firecloud_api(
            firecloud_api.list_repository_methods,
            namespace=self.repo_namespace,
            name=self.repo_method_name,
        )

        snapshots.sort(key=lambda x: x["snapshotId"], reverse=True)
        return snapshots

    def delete_old_method_snapshots(self, nkeep: int) -> None:
        """
        Delete all but `n` of the most recent snapshots of the method. This might fail
        if the service account doesn't have OWNER permission on the method configuration
        namespace. This can't be set using the firecloud package, so POST to the
        `setConfigNamespaceACL` endpoint on https://api.firecloud.org/ once to resolve
        this problem.

        :param nkeep: the number of snapshots to keep
        """

        snapshots = self.get_method_snapshots()

        to_delete = snapshots[nkeep:]
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
        self, entity_type: str, pandera_schema: Type[PanderaBaseSchema]
    ) -> TypedDataFrame[PanderaBaseSchema]:
        """
        Get a data frame of entities from a Terra data table.

        :param entity_type: the kind of entity (e.g. "sample")
        :param pandera_schema: a Pandera schema for the output data frame
        """

        echo(f"Getting {entity_type} entities")
        j = call_firecloud_api(
            firecloud_api.get_entities,
            namespace=self.workspace_namespace,
            workspace=self.workspace_name,
            etype=entity_type,
        )

        records = [{f"{entity_type}_id": x["name"], **x["attributes"]} for x in j]

        return TypedDataFrame[pandera_schema](pd.DataFrame(records))

    def upload_entities(self, df: pd.DataFrame) -> None:
        """
        Upload a data frame of entities to a Terra data table.

        :param df: a data frame of entities
        """

        echo(f"{len(df)} entities to upload to Terra")

        for batch in batch_evenly(df, max_batch_size=500):
            with tempfile.NamedTemporaryFile(suffix="tsv") as f:
                batch.to_csv(f, sep="\t", index=False)  # pyright: ignore
                f.flush()

                echo(f"Upserting {len(batch)} entities to Terra")
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

        echo("Creating workspace config")
        call_firecloud_api(
            firecloud_api.create_workspace_config,
            namespace=self.workspace_namespace,
            workspace=self.workspace_name,
            body=config_body,
        )

        echo("Setting workspace config ACL")
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

        echo("Update workspace config")
        call_firecloud_api(
            firecloud_api.update_workspace_config,
            namespace=self.workspace_namespace,
            workspace=self.workspace_name,
            cnamespace=terra_workflow.repo_namespace,
            configname=terra_workflow.repo_method_name,
            body=config_body,
        )

        echo("Setting workspace config ACL")
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

        terra_workflow.persist_method_on_gcs()
        assert terra_workflow.persisted_wdl_script is not None

        echo("Creating method")
        snapshot = call_firecloud_api(
            firecloud_api.update_repository_method,
            namespace=terra_workflow.repo_namespace,
            method=terra_workflow.repo_method_name,
            synopsis=terra_workflow.method_synopsis,
            wdl=terra_workflow.workflow_wdl_path,
        )

        echo("Setting method ACL")
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
        assert terra_workflow.persisted_wdl_script is not None

        with tempfile.NamedTemporaryFile("w") as f:
            f.write(terra_workflow.persisted_wdl_script["wdl"])
            f.flush()

            echo("Updating method")
            snapshot = call_firecloud_api(
                firecloud_api.update_repository_method,
                namespace=terra_workflow.repo_namespace,
                method=terra_workflow.repo_method_name,
                synopsis=terra_workflow.method_synopsis,
                wdl=f.name,
            )

        echo("Setting method ACL")
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

        # inject the workflow version and URL into inputs so it gets stored in job
        # submissions
        assert terra_workflow.persisted_wdl_script is not None

        if "version" in terra_workflow.persisted_wdl_script:
            terra_workflow.method_config["inputs"][
                f"{terra_workflow.repo_method_name}.workflow_version"
            ] = f'"{terra_workflow.persisted_wdl_script["version"]}"'

        terra_workflow.method_config["inputs"][
            f"{terra_workflow.repo_method_name}.workflow_url"
        ] = f'"{terra_workflow.persisted_wdl_script["public_url"]}"'

        echo("Checking for existing workspace config")
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
        terra_workflow.delete_old_method_snapshots(nkeep=20)

    def submit_workflow_run(
        self, terra_workflow: TerraWorkflow, **kwargs: Unpack[TerraJobSubmissionKwargs]
    ) -> None:
        """
        Submit a run of a workflow.

        :param terra_workflow: a `TerraWorkflow` instance
        """

        echo(f"Submitting {terra_workflow.repo_method_name} job")
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
                .rstrip("+00:00")
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

    def collect_workflow_outputs(
        self, since: datetime.datetime | None = None
    ) -> list[task_result_insert_input]:
        """
        Collect workflow outputs and metadata about the jobs that ran them.

        :param since: don't collect outputs for job submissions before this `datetime`
        :return: a list of Gumbo `task_result` objects to insert
        """

        echo("Getting previous job submissions")
        submissions = pd.DataFrame(
            call_firecloud_api(
                firecloud_api.list_submissions,
                namespace=self.workspace_namespace,
                workspace=self.workspace_name,
            )
        ).convert_dtypes()

        submissions = expand_dict_columns(submissions).convert_dtypes()

        if since is not None:
            # the list of submissions can grow very long and the endpoint doesn't
            # support filtering
            submissions = submissions.loc[submissions["submissionDate"].ge(since)]

        # only collect outputs from submissions with successful workflow runs
        submissions = submissions.loc[submissions["workflowStatuses__Succeeded"].gt(0)]

        outputs = []

        for _, s in submissions.iterrows():
            # start constructing a common output object to be used for `task_result`
            # inserts
            base_o = task_result_insert_input(
                terra_method_config_name=str(s["methodConfigurationName"]),
                terra_method_config_namespace=str(s["methodConfigurationNamespace"]),
                terra_submission_id=str(s["submissionId"]),
                terra_workspace_name=self.workspace_name,
                terra_workspace_namespace=self.workspace_namespace,
            )

            # get the workflows for this job submission (often there is only 1)
            echo(f"Getting workflows for job submission {s["submissionId"]}")
            submission = call_firecloud_api(
                firecloud_api.get_submission,
                namespace=self.workspace_namespace,
                workspace=self.workspace_name,
                submission_id=s["submissionId"],
            )

            for w in submission["workflows"]:
                if "workflowId" not in w:
                    # this workflow didn't manage to start
                    continue

                base_o.terra_workflow_id = w["workflowId"]
                base_o.terra_entity_name = w["workflowEntity"]["entityName"]
                base_o.terra_entity_type = w["workflowEntity"]["entityType"]
                base_o.created_at = pd.Timestamp(
                    ts_input=w["statusLastChangedDate"]
                ).isoformat()  # pyright: ignore

                echo(f"Getting workflow {w["workflowId"]} metadata")
                wmd = call_firecloud_api(
                    firecloud_api.get_workflow_metadata,
                    namespace=self.workspace_namespace,
                    workspace=self.workspace_name,
                    submission_id=s["submissionId"],
                    workflow_id=w["workflowId"],
                    include_key=[
                        "inputs",
                        "labels",
                        "outputs",
                        "status",
                        "workflowName",
                        "workflowRoot",
                    ],
                )

                if wmd["status"] != "Succeeded":
                    continue

                base_o.workflow_name = wmd["workflowName"]
                base_o.terra_workflow_inputs = wmd["inputs"]
                base_o.terra_workflow_root_dir = wmd["workflowRoot"]
                base_o.terra_workspace_id = wmd["labels"]["workspace-id"]

                # workflow URL and version should be injected into the inputs when the
                # method config is deployed
                if "workflow_source_url" in wmd["inputs"]:
                    base_o.workflow_source_url = wmd["inputs"]["workflow_source_url"]

                if "workflow_version" in wmd["inputs"]:
                    base_o.workflow_version = wmd["inputs"]["workflow_version"]

                # iterate through the outputs and make a `task_result` object for each
                for label, output in wmd["outputs"].items():
                    # all attributes up to now have been invariant between outputs
                    o = base_o.copy()

                    # the workflow outputs are named like `workflow_name.output_name`
                    o.label = label.rsplit(".", maxsplit=1)[-1]

                    if isinstance(output, list):
                        # we shouldn't be writing any workflows that output lists, so
                        # these would just be intermediate or legacy outputs
                        continue
                    elif isinstance(output, str):
                        if output.startswith("gs://"):
                            o.url = output
                            o.format = pathlib.Path(o.url).suffix[1:].upper()
                        else:
                            # it's a plain string
                            o.value = {"value": output}
                    else:
                        # it's a number, bool, or something else
                        o.value = {"value": output}

                    outputs.append(o)

        return outputs


def call_firecloud_api(func: Callable, *args: Any, **kwargs: Any) -> Any:
    """
    Call a Firecloud API endpoint and check the response for a valid HTTP status code.

    :param func: a `firecloud.api` method
    :param args: arguments to `func`
    :param kwargs: keyword arguments to `func`
    :return: the API response, if any
    """

    res = maybe_retry(
        func,
        retryable_exceptions=(requests.ConnectionError, requests.ConnectTimeout),
        max_retries=4,
        *args,
        **kwargs,
    )

    if 200 <= res.status_code <= 299:
        try:
            return res.json()
        except requests.JSONDecodeError:
            return res.text

    try:
        raise requests.RequestException(f"HTTP {res.status_code} error: {res.json()}")
    except Exception as e:
        # it's returning HTML or we can't parse the JSON
        echo(f"Error getting response as JSON: {e}", err=True)
        echo(f"Response text: {res.text}", err=True)
        raise requests.RequestException(f"HTTP {res.status_code} error")
