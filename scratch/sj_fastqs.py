import gzip
import json
from typing import Any

import gcsfs
import pandas as pd
from google.cloud import storage
from nebelung.terra_workspace import TerraWorkspace
from pd_flatten import pd_flatten

pd.set_option("display.max_columns", 30)
pd.set_option("display.max_colwidth", 50)
pd.set_option("display.max_info_columns", 30)
pd.set_option("display.max_info_rows", 20)
pd.set_option("display.max_rows", 20)
pd.set_option("display.max_seq_items", None)
pd.set_option("display.width", 200)
pd.set_option("expand_frame_repr", True)
pd.set_option("mode.chained_assignment", "warn")


def list_blobs(
    bucket_name: str, prefix: str | None = None, glob: str | None = None
) -> pd.DataFrame:
    """
    Get the names and sizes of existing blobs in a GCS bucket.

    :param bucket_name: the name of the GCS bucket
    :param prefix: an optional prefix for listing
    :param glob: an optional glob for listing
    :return: a data frame of object names and sizes
    """

    storage_client = storage.Client(project="depmap-omics")

    if prefix is not None and glob is not None:
        raise ValueError("At most one of `prefix` and `glob` can be specified")
    elif prefix is not None:
        pages = storage_client.list_blobs(
            bucket_or_name=bucket_name,
            prefix=prefix,
            delimiter="/",
            fields="items(name,crc32c,size,updated),nextPageToken",
        ).pages
    elif glob is not None:
        pages = storage_client.list_blobs(
            bucket_or_name=bucket_name,
            match_glob=glob,
            fields="items(name,crc32c,size,updated),nextPageToken",
        ).pages
    else:
        pages = storage_client.list_blobs(
            bucket_or_name=bucket_name,
            fields="items(name,crc32c,size,updated),nextPageToken",
        ).pages

    blobs = []

    for page in pages:
        blobs.extend(
            [
                {
                    "url": "gs://" + bucket_name + "/" + x.name,
                    "crc32c": x.crc32c,
                    "size": x.size,
                    "gcs_obj_updated_at": x.updated,
                }
                for x in page
            ]
        )

    return pd.DataFrame(blobs)


def make_readgroup_header(sample_id: str, fastq_url: str) -> str:
    """
    Read the first header line from `fastq_url` (a .fastq.gz in GCS)
    and construct an @RG header suitable for BWA.
    """

    fs = gcsfs.GCSFileSystem()

    # Read the very first line of the FASTQ
    with fs.open(fastq_url, "rb") as f:
        with gzip.GzipFile(fileobj=f) as gz:
            header = gz.readline().decode("utf-8").strip()

    # Example header:
    # @LH00297:152:22WKHHLT4:7:1101:9875:1042 1:N:0:TAATTTCACT+TGAGTCGGTT
    if not header.startswith("@"):
        raise ValueError(f"Unexpected FASTQ header line: {header}")

    # Split into left (instrument..coords) and right (read/filter/index info)
    print(header)
    left, right = header[1:].split(" ", 1)
    parts = left.split(":")
    if len(parts) < 4:
        raise ValueError(f"Malformed Illumina header: {header}")

    flowcell = parts[2]
    lane = parts[3]

    # Parse index from right side (READ:FILT:CTRL:INDEX)
    r_parts = right.split(":")
    if len(r_parts) < 4:
        raise ValueError(f"Malformed Illumina header index section: {header}")
    index_seq = r_parts[3]

    # Build RG fields
    rg_id = f"{flowcell}.L{lane}"
    rg_pu = f"{flowcell}.L{lane}.{index_seq}"

    # Construct BWA-ready @RG line (tabs required)
    rg = f"@RG\tID:{rg_id}\tSM:{sample_id}\tPU:{rg_pu}\tPL:ILLUMINA"
    return rg


df = list_blobs(bucket_name="depmap-stjude")
assert bool(~df["size"].duplicated().any())
df = df.loc[:, ["url"]]

df["name"] = df["url"].str.split("/").str.get(-1)
df["sample_id"] = df["name"].str.extract(r"^\d+_(.+)_S\d+_L\d+_R\d+_\d+.fastq.gz$")

df["rg"] = df.apply(lambda x: make_readgroup_header(x["sample_id"], x["url"]), axis=1)

df.to_parquet("./data/sj_fastqs.parquet", index=False)
df = pd.read_parquet("./data/sj_fastqs.parquet")

df["strand"] = "forward"
df.loc[df["name"].str.contains("_R2_"), "strand"] = "reverse"


def make_sample(g: pd.DataFrame) -> dict[str, Any]:
    return {
        "readgroup_id": g["sample_id"].iloc[0],
        "readgroup_name": g["rg"].iloc[0],
        "forward_fastq": g.loc[g["strand"].eq("forward"), "url"].squeeze(),
        "reverse_fastq": g.loc[g["strand"].eq("reverse"), "url"].squeeze(),
    }


samples = df.copy()

samples = samples.merge(
    samples.groupby(["sample_id", "rg"])
    .apply(make_sample)
    .to_frame(name="fastqs_pe")
    .reset_index(),
    how="inner",
    on=["sample_id", "rg"],
)[["sample_id", "fastqs_pe"]]

# samples = (
#     samples.groupby("sample_id")["fastqs_pe"]
#     .agg(lambda x: json.dumps(list(x)))
#     .reset_index()
# )

samples = pd_flatten(samples)

samples = (
    samples.groupby("sample_id")[
        [
            "fastqs_pe__readgroup_id",
            "fastqs_pe__readgroup_name",
            "fastqs_pe__forward_fastq",
            "fastqs_pe__reverse_fastq",
        ]
    ]
    .agg(lambda x: json.dumps(list(x)))
    .reset_index()
)

# wgs_samples = TerraWorkspace("broad-firecloud-ccle", "depmap-omics-wgs").get_entities(
#     "sample"
# )

ref_vals = (
    wgs_samples.loc[
        :,
        [
            "ref_alt",
            "ref_amb",
            "ref_ann",
            "ref_bwt",
            "ref_dict",
            "ref_fasta",
            "ref_fasta_index",
            "ref_pac",
            "ref_sa",
        ],
    ]
    .iloc[0]
    .to_dict()
)

samples = samples.assign(**ref_vals)

sample_ids = samples.pop("sample_id")
samples.insert(0, "entity:sample_id", sample_ids)

samples = samples.rename(
    columns={
        "fastqs_pe__readgroup_id": "fastqs_pe_readgroup_id",
        "fastqs_pe__readgroup_name": "fastqs_pe_readgroup_name",
        "fastqs_pe__forward_fastq": "fastqs_pe_forward",
        "fastqs_pe__reverse_fastq": "fastqs_pe_reverse",
    }
)

samples["fastq_se"] = pd.NA
samples["fastq_se_readgroup"] = pd.NA
samples["fastq_se_readgroup_id"] = pd.NA
samples["input_type"] = "FASTQ"

TerraWorkspace("broad-firecloud-ccle", "depmap-sj-wgs-align").upload_entities(samples)
