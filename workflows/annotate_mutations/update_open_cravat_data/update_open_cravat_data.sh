#!/bin/bash

set -euo pipefail

oc module install-base
oc module install -y \
  vcfreporter \
  brca1_func_assay \
  ccre_screen \
  gtex \
  gwas_catalog \
  pharmgkb \
  provean \
  revel \
  spliceai

RSYNC_SRC="$(oc config md)"
RSYNC_DEST="${GS_URL_PREFIX}/$(date -I)"

echo "Syncing ${RSYNC_SRC} to ${RSYNC_DEST}"
gcloud config set storage/parallel_composite_upload_enabled True
gcloud storage rsync \
    --recursive \
    --delete-unmatched-destination-objects \
    --no-user-output-enabled \
    "${RSYNC_SRC}" \
    "${RSYNC_DEST}"

echo "Done."
