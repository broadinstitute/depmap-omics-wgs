import base64
import json
import logging

import functions_framework
import google.cloud.logging
from cloudevents.http import CloudEvent

from omics_wgs_pipeline.terra import refresh_terra_samples


@functions_framework.cloud_event
def run(cloud_event: CloudEvent) -> None:
    """
    Bypass omics_wgs_pipeline CLI commands to call the underlying functions directly (needed for
    remote execution inside a GCP Function).

    :param cloud_event: the pub/sub CloudEvent payload
    """

    client = google.cloud.logging.Client()
    client.setup_logging(log_level=logging.INFO)

    ce_data = json.loads(base64.b64decode(cloud_event.data["message"]["data"]).decode())

    refresh_terra_samples(
        repo_namespace=ce_data["repo_namespace"],
        repo_method_name_preprocess=ce_data["repo_method_name_preprocess"],
        workspace_namespace=ce_data["workspace_namespace"],
        workspace_name=ce_data["workspace_name"],
        gumbo_env=ce_data["gumbo_env"],
        use_default_service_account=True,
    )
