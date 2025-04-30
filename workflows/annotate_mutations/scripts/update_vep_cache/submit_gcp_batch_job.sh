#!/bin/zsh

set -euo pipefail

gcloud storage cp \
  ./update_vep.sh \
  gs://ccleparams/scripts/update_vep.sh

gcloud batch jobs submit \
    --job-prefix update-vep \
    --config=./gcp_batch_job_config.json \
    --location=us-central1
