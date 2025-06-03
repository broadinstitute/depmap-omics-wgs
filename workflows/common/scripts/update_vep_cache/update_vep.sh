#!/bin/bash

set -euo pipefail

curl -O https://ftp.ensembl.org/pub/release-113/variation/indexed_vep_cache/homo_sapiens_vep_113_GRCh38.tar.gz
tar xzf homo_sapiens_vep_113_GRCh38.tar.gz

BASE_DIR="homo_sapiens/113_GRCh38"
RSYNC_SRC="/tmp/vep_out"

mkdir -p "$RSYNC_SRC"

# Loop through all subdirectories in BASE_DIR
for dir in "$BASE_DIR"/*/; do
    subdir=$(basename "$dir")
    archive_path="${RSYNC_SRC}/${subdir}.gz"

    # Create archive with:
    # - all top-level files in BASE_DIR
    # - the current subdirectory
    (
        cd homo_sapiens || exit 1

        find "113_GRCh38" -maxdepth 1 -type f -printf "%p\n" | \
          xargs tar -czf "$archive_path" "113_GRCh38/$subdir"
    )

    echo "Created archive: $archive_path"
done

RSYNC_DEST="${GS_URL_PREFIX}/v113"

echo "Syncing ${RSYNC_SRC} to ${RSYNC_DEST}"
gcloud config set storage/parallel_composite_upload_enabled True
gcloud storage rsync \
    --recursive \
    --delete-unmatched-destination-objects \
    --no-user-output-enabled \
    "${RSYNC_SRC}" \
    "${RSYNC_DEST}"

echo "Done."
