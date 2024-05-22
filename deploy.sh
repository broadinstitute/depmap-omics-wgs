#!/bin/zsh

set -euo pipefail

pre-commit run poetry-export --all-files

echo "Updating Terra workflow (preprocess)"
poetry run python -m omics_wgs_pipeline update-workflow \
  --config-path config.toml \
  update-workflow \
  --workflow-name preprocess_wgs_sample

echo "Updating GCP Function"
gcloud functions deploy omics_wgs_pipeline \
  --gen2 \
  --runtime="python312" \
  --region="us-central1" \
  --source=. \
  --run-service-account="omics-pipeline-runner@depmap-omics.iam.gserviceaccount.com" \
  --entry-point="run" \
  --trigger-topic="run-omics-pipeline" \
  --timeout=600 \
  --memory="4GB" \
  --cpu=1 \
  --max-instances=1 \
  --concurrency=1
