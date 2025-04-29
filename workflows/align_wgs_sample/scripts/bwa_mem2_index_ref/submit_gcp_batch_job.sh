#!/bin/zsh

set -euo pipefail

gcloud storage cp \
  ./bwa_mem2_index_ref.sh \
  gs://cds-pipelines/scripts/bwa_mem2_index_ref.sh

gcloud batch jobs submit \
    --job-prefix bwa-mem2-index-ref \
    --config=./gcp_batch_job_config.json \
    --location=us-central1
