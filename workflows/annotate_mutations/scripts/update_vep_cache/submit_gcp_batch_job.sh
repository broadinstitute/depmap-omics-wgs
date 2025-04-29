#!/bin/zsh

set -euo pipefail

gcloud storage cp \
  ./update_open_cravat_data.sh \
  gs://cds-pipelines/scripts/update_open_cravat_data.sh

gcloud batch jobs submit \
    --job-prefix update-open-cravat-data \
    --config=./gcp_batch_job_config.json \
    --location=us-central1
