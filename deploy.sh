#!/bin/zsh

set -euo pipefail

pre-commit run poetry-export --all-files

GCR_LOCATION="us-docker.pkg.dev/depmap-omics/public/omics_wgs_pipeline"
PKG_VERSION="0.1.0"
IMAGE_TAG="${GCR_LOCATION}:${PKG_VERSION}"

echo "Building Docker image"
docker build . -t "${IMAGE_TAG}" --platform="linux/amd64"
echo "Pushing Docker image to GCR"
gcloud auth configure-docker us-docker.pkg.dev
docker push "${IMAGE_TAG}"

echo "Deleting old GCR images"
OLD_DIGESTS=$(gcloud container images list-tags ${GCR_LOCATION} \
  --sort-by="~timestamp" --format=json | jq -r '.[2:] | .[].digest')

if [ -n "${OLD_DIGESTS}" ]; then
  while read -r DIGEST; do
    gcloud container images delete "${GCR_LOCATION}@${DIGEST}" --quiet
  done <<< "${OLD_DIGESTS}"
fi

echo "Getting latest image hash"
LATEST_DIGEST=$(gcloud container images list-tags ${GCR_LOCATION} \
  --sort-by="~timestamp" --limit=1 --format=json | jq -r  '.[].digest')

echo "Updating production tag"
gcloud container images add-tag --quiet \
  "${GCR_LOCATION}@${LATEST_DIGEST}" \
  "${GCR_LOCATION}:production"

echo "Updating Terra workflow (preprocess)"
poetry run python -m omics_wgs_pipeline update-workflow \
  --repo-namespace="cds-pipelines" \
  --repo-method-name="preprocess" \
  --workspace-namespace="broad-firecloud-ccle" \
  --workspace-name="omics_wgs_pipeline" \
  --method-config-name="preprocess" \
  --method-synopsis="Align and preprocess a single sample" \
  --workflow-wdl-path="./workflows/preprocess/preprocess.wdl" \
  --method-config-json-path="./workflows/preprocess/preprocess.json" \
  --firecloud-owners="dmccabe@broadinstitute.org" \
  --firecloud-owners="wgs-omics-pipeline-runner@depmap-omics.iam.gserviceaccount.com" \
  --image-hash="${LATEST_DIGEST}"

echo "Updating GCP Function"
gcloud functions deploy omics_wgs_pipeline \
  --gen2 \
  --runtime="python39" \
  --region="us-central1" \
  --source=. \
  --run-service-account="wgs-omics-pipeline-runner@depmap-omics.iam.gserviceaccount.com" \
  --entry-point="run" \
  --trigger-topic="run-wgs-omics-pipeline" \
  --timeout=600 \
  --memory="8GB" \
  --cpu=2 \
  --max-instances=1 \
  --concurrency=1
