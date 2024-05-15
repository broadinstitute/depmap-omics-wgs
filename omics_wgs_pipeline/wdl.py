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
    """
    Save a copy of a local WDL script and its dependencies to GCS, rewriting the import
    statements to refer to the public GCS URLs.

    :param bucket: a GCS bucket
    :param wdl_path: the absolute path to a WDL script
    :param subpath: a prefix to prepend to the destiation URL
    :return: a dictionary of the uploaded WDL and its public GCS URL
    """

    with open(wdl_path, "r") as f:
        wdl_lines = f.readlines()

    wdl_basename = os.path.basename(wdl_path)
    versioned_subpath = subpath
    buffer = []

    # include the version number on the uploaded object path if one is found
    for line in wdl_lines:
        if version_match := re.match(VERSION_PATTERN, line):
            versioned_subpath = "/".join(
                [versioned_subpath, Path(wdl_basename).stem, version_match[1]]
            )

    for line in wdl_lines:
        if import_match := re.match(IMPORT_PATTERN, line):
            # need to upload the dependent WDL file to GCS and rewrite the `import`
            # statement in this WDL file to the dependent file's public GCS URL

            rel_path = import_match[1]  # relative path to the dependent file
            import_alias = import_match[2]  # the `as <x>` component of the statement

            # absolute path to the dependent file
            abs_wdl_path = Path.joinpath(wdl_path.parent.absolute(), rel_path)

            # append the local relative path to the end of parent's prefix
            imported_versioned_subpath = "/".join(
                [versioned_subpath, Path(rel_path).parent.name]
            )

            # recurse: upload the dependent WDL and get its URL
            imported_public_gcs_url = persist_wdl_script(
                bucket, abs_wdl_path, imported_versioned_subpath
            )["public_url"]

            # replace the local relative import with the absolute one
            converted_line = f'import "{imported_public_gcs_url}" as {import_alias}'
            buffer.append(converted_line)

        else:
            buffer.append(line.rstrip())

    # construct the final version of this WDL script
    converted_wdl = "\n".join(buffer)

    # upload the script and make it public
    wdl_gs_path = "/".join([versioned_subpath, wdl_basename])
    blob = bucket.blob(wdl_gs_path)
    blob.upload_from_string(converted_wdl)
    blob.make_public()

    return {"wdl": converted_wdl, "public_url": blob.public_url}
