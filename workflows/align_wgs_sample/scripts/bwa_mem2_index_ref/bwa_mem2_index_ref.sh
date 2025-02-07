#!/bin/bash

set -euo pipefail

mkdir -p out && cd out
LOCAL_REF_FASTA=$(basename "${REF_FASTA}")
gcloud storage cp "${REF_FASTA}" "${LOCAL_REF_FASTA}"

bwa-mem2 index "${LOCAL_REF_FASTA}"

RSYNC_SRC="."
RSYNC_DEST="${GCS_FOLDER}"

echo "Syncing ${RSYNC_SRC} to ${RSYNC_DEST}"
gcloud config set storage/parallel_composite_upload_enabled True
gcloud storage rsync \
    --delete-unmatched-destination-objects \
    "${RSYNC_SRC}" \
    "${RSYNC_DEST}"

echo "Done."
