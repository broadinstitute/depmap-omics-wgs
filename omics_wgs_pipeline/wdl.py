import os
import re
from pathlib import Path

from google.cloud import storage

from omics_wgs_pipeline.types import PersistedWdl

IMPORT_PATTERN = re.compile(r"^import\s+\"(?!http)([^\"]+)\"\s+as\s+(\S+)$")
VERSION_PATTERN = re.compile(r"^version\s+([\d.]+)$")


def persist_wdl_script(
    bucket: storage.Bucket, wdl_path: Path, subpath: str = ""
) -> PersistedWdl:
    with open(wdl_path, "r") as f:
        wdl_lines = f.readlines()

    wdl_basename = os.path.basename(wdl_path)
    versioned_subpath = subpath
    buffer = []

    for line in wdl_lines:
        if version_match := re.match(VERSION_PATTERN, line):
            versioned_subpath = "/".join(
                [versioned_subpath, Path(wdl_basename).stem, version_match[1]]
            )

    for line in wdl_lines:
        if import_match := re.match(IMPORT_PATTERN, line):
            rel_path = import_match[1]
            import_alias = import_match[2]

            abs_wdl_path = Path.joinpath(wdl_path.parent.absolute(), rel_path)
            imported_versioned_subpath = "/".join(
                [versioned_subpath, Path(rel_path).parent.name]
            )
            imported_public_gcs_url = persist_wdl_script(
                bucket, abs_wdl_path, imported_versioned_subpath
            )["public_url"]

            converted_line = f'import "{imported_public_gcs_url}" as {import_alias}'
            buffer.append(converted_line)

        else:
            buffer.append(line.rstrip())

    converted_wdl = "\n".join(buffer)

    wdl_gs_path = "/".join([versioned_subpath, wdl_basename])
    blob = bucket.blob(wdl_gs_path)
    blob.upload_from_string(converted_wdl)
    blob.make_public()

    return {"wdl": converted_wdl, "public_url": blob.public_url}
