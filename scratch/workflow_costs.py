import os
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from pprint import pp

import pandas as pd
from firecloud_api_cds import api as firecloud_api
from nebelung.terra_workspace import TerraWorkspace
from nebelung.utils import call_firecloud_api
from pd_flatten import pd_flatten


class BoundedThreadPoolExecutor:
    def __init__(
        self,
        max_workers: int | None = None,
        max_queue: int | None = None,
    ):
        """
        Use a semaphore to bound the number of in-flight tasks in a thread pool.

        :param max_workers: max number of workers to use
        :param max_queue: number of jobs to keep in the task queue
        """

        if max_workers is None:
            max_workers = min(32, (os.cpu_count() or 1) + 4)
        if max_workers <= 0:
            raise ValueError("max_workers must be greater than 0")

        if max_queue is None:
            max_queue = 4 * max_workers  # good for IO-bound, otherwise use 1:1 ratio
        if max_queue <= 0:
            raise ValueError("max_queue must be greater than 0")

        self._executor = ThreadPoolExecutor(max_workers=max_workers)
        self._sem = threading.Semaphore(max_queue)

    def submit(self, fn, *args, **kwargs):
        # Apply backpressure
        self._sem.acquire()
        fut = self._executor.submit(fn, *args, **kwargs)
        fut.add_done_callback(lambda _: self._sem.release())
        return fut

    def shutdown(self, wait: bool = True):
        self._executor.shutdown(wait=wait)

    # --- context manager support ---
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.shutdown(wait=True)
        return False


# workspace_namespace = "broad-firecloud-ccle"
# workspace_name = "depmap-omics-wgs"
# date_start = "2025-09-01"
# date_end = "2025-11-30"

workspace_namespace = "broad-firecloud-ccle"
workspace_name = "DepMap_WGS_CN"
date_start = "2024-07-03"
date_end = "2024-11-05"

submissions = pd.DataFrame(
    call_firecloud_api(
        firecloud_api.list_submissions,
        namespace=workspace_namespace,
        workspace=workspace_name,
    )
)

submissions = pd_flatten(submissions).convert_dtypes()

submissions["submissionDate"] = pd.to_datetime(submissions["submissionDate"])

submissions = submissions.loc[
    submissions["submissionDate"].ge(date_start)
    & submissions["submissionDate"].lt(date_end)
]

submissions = submissions.loc[
    (
        submissions["methodConfigurationName"].str.startswith("WGS_preprocessing")
        | submissions["methodConfigurationName"].str.startswith("WGS_pipeline")
        | submissions["methodConfigurationName"].str.startswith("guide_mutation_binary")
        | submissions["methodConfigurationName"].str.startswith("vcf_to_depmap")
        | submissions["methodConfigurationName"].str.startswith("vep_sv")
        | submissions["methodConfigurationName"].str.startswith("Manta_SomaticSV_v1_0")
        | submissions["methodConfigurationName"].str.startswith("PureCN_dev")
    )
    & ~submissions["submissionEntity__entityName"].isin(
        ["24Q4_DRAGEN_GIAB", "24Q4_DRAGEN_GIAB", "24Q4_NON_DRAGEN_GIAB"]
    )
]

# submissions = submissions.sample(n=10)


def get_submission_costs(s):
    sid = s["submissionId"]

    # start constructing a common output object to be used for `task_result`
    # inserts
    base_wf_data = dict(
        terra_submission_id=str(sid),
        terra_workspace_name=workspace_name,
        terra_workspace_namespace=workspace_namespace,
        method_configuration_name=s["methodConfigurationName"],
    )

    # get the workflows for this job submission (often there is only 1)
    print(f"Getting workflows for job submission {sid}")
    submission = call_firecloud_api(
        firecloud_api.get_submission,
        namespace=workspace_namespace,
        workspace=workspace_name,
        submission_id=sid,
    )

    arr = []

    for w in submission["workflows"]:
        if "workflowId" not in w:
            # this workflow didn't manage to start
            continue

        if "workflowEntity" not in w:
            # this workflow didn't manage to start
            continue

        wf_data = base_wf_data.copy()

        wid = w["workflowId"]
        wf_data["cost"] = w["cost"] if "cost" in w else 0
        wf_data["status"] = w["status"]
        wf_data["terra_workflow_id"] = wid
        wf_data["terra_entity_name"] = w["workflowEntity"]["entityName"]
        wf_data["terra_entity_type"] = w["workflowEntity"]["entityType"]

        if workspace_name == "depmap-omics-wgs":
            wmd = call_firecloud_api(
                firecloud_api.get_workflow_metadata,
                namespace=workspace_namespace,
                workspace=workspace_name,
                submission_id=sid,
                workflow_id=wid,
                include_key=["workflowName"],
            )

        wf_data["workflow_name"] = wmd["workflowName"]

        arr.append(wf_data)

    return arr


dfs = []

with BoundedThreadPoolExecutor(max_workers=10, max_queue=100) as executor:
    futures = {
        executor.submit(get_submission_costs, s) for _, s in submissions.iterrows()
    }

    for future in as_completed(futures):
        dfs.extend(future.result())

df = pd.DataFrame.from_records(dfs)

if workspace_name == "DepMap_WGS_CN":
    ws = TerraWorkspace("broad-firecloud-ccle", "DepMap_WGS_CN")
    sample_sets = ws.get_entities("sample_set")
    sample_ids = [
        x["entityName"]
        for x in sample_sets.loc[
            sample_sets["sample_set_id"].eq("24Q4"), "samples"
        ].squeeze()["items"]
    ]

    df = df.loc[df["terra_entity_name"].isin(sample_ids)]

df.to_parquet(f"./data/costs/{workspace_name}.parquet")
