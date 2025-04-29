import base64
import json
import logging

import functions_framework
import google.cloud.logging
from cloudevents.http import CloudEvent


@functions_framework.cloud_event
def run(cloud_event: CloudEvent) -> None:
    """
    Bypass depmap_omics_wgs CLI commands to call the underlying functions directly (needed for
    remote execution inside a GCP Function).

    :param cloud_event: the pub/sub CloudEvent payload
    """

    client = google.cloud.logging.Client()
    client.setup_logging(log_level=logging.INFO)

    ce_data = json.loads(base64.b64decode(cloud_event.data["message"]["data"]).decode())

    raise NotImplementedError
