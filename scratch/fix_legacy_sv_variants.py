import itertools
import tomllib
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Iterable

import pandas as pd
from firecloud_api_cds import api as firecloud_api
from google.cloud import storage
from nebelung.terra_workspace import TerraWorkspace
from nebelung.utils import call_firecloud_api

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

    storage_client = storage.Client()

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


def filter_svs(
    df: pd.DataFrame,
    cosmic_fusion_gene_pairs_path: Path,
    onco_tsg_path: Path,
    min_depth: int,
    sv_gnomad_cutoff: float,
    large_sv_size: int,
) -> pd.DataFrame:
    """
    Filter structural variants while rescuing clinically important ones.

    Applies size and frequency filters but rescues large SVs, variants affecting
    oncogenes/tumor suppressors, and variants creating known COSMIC fusion pairs.

    :param df: SVs in BEDPE format
    :param cosmic_fusion_gene_pairs_path: Path to file containing COSMIC fusion gene
        pairs
    :param onco_tsg_path: Path to file containing oncogene and tumor suppressor gene
        symbols
    :param min_depth: Minimum read depth (PR+SR) of variant
    :param sv_gnomad_cutoff: Maximum gnomAD allele frequency for an SV to be considered
        somatic
    :param large_sv_size: Size threshold beyond which SVs are considered large and need
        to be rescued
    :returns: Filtered SVs in BEDPE format
    """

    df["SVLEN_A"] = df["SVLEN_A"].astype("Int64")

    # drop variants below minimum read depth
    df["dp"] = df["PR"].str.split(",", expand=True).astype("Int64").sum(axis=1) + df[
        "SR"
    ].str.split(",", expand=True).astype("Int64").sum(axis=1)
    df = df.loc[df["dp"].ge(min_depth)]

    # drop variants shorter than 50
    df = df.loc[df["SVLEN_A"].isna() | df["SVLEN_A"].abs().ge(50)].copy()

    onco_tsg_df = pd.read_csv(
        onco_tsg_path, sep="\t", dtype="string", usecols=["hugo_symbol"]
    )
    oncogenes_and_ts = set(onco_tsg_df["hugo_symbol"])

    cosmic = pd.read_csv(cosmic_fusion_gene_pairs_path, dtype="string")
    cosmic_pairs = list(
        zip(
            list(cosmic["five_prime_gene_symbol"]),
            list(cosmic["three_prime_gene_symbol"]),
        )
    )
    cosmic_pairs_sorted = set([tuple(sorted(elem)) for elem in cosmic_pairs])

    df["Rescue"] = False

    # rescue large SVs
    df.loc[df["SVLEN_A"].abs().ge(large_sv_size), "Rescue"] = True

    # rescue breakpoints that fall on oncogenes or tumor suppressors
    df["onco_ts_overlap_A"] = df["SYMBOL_A"].apply(
        onco_ts_overlap, oncogenes_and_ts=oncogenes_and_ts
    )
    df["onco_ts_overlap_B"] = df["SYMBOL_B"].apply(
        onco_ts_overlap, oncogenes_and_ts=oncogenes_and_ts
    )
    df.loc[df["onco_ts_overlap_A"] | df["onco_ts_overlap_B"], "Rescue"] = True

    # rescue gene pairs in cosmic
    df["pair_in_cosmic"] = df.apply(
        lambda row: list_all_pairs(
            row["SYMBOL_A"],
            row["SYMBOL_B"],
            cosmic_pairs_sorted,
        ),
        axis=1,
    )
    df.loc[df["pair_in_cosmic"], "Rescue"] = True

    # gnomad AF parsing
    df["max_af"] = (
        df["vep_SV_overlap_AF_A"]
        .fillna("")
        .str.split("&")
        .apply(lambda x: max([float(e) if e != "" else 0 for e in x]))
    )

    # filter while keeping rescues
    df = df.loc[
        df["Rescue"] | df["max_af"].lt(sv_gnomad_cutoff),
        [
            "CHROM_A",
            "START_A",
            "END_A",
            "ID",
            "STRAND_A",
            "TYPE",
            "FILTER",
            "REF_A",
            "ALT_A",
            "SVLEN_A",
            "MATEID_A",
            "SVINSLEN_A",
            "BND_DEPTH_A",
            "MATE_BND_DEPTH_A",
            "SYMBOL_A",
            "GENEID_A",
            "vep_SV_overlap_name_A",
            "vep_SV_overlap_AF_A",
            "CHROM_B",
            "START_B",
            "END_B",
            "STRAND_B",
            "REF_B",
            "ALT_B",
            "SVLEN_B",
            "MATEID_B",
            "SVINSLEN_B",
            "BND_DEPTH_B",
            "MATE_BND_DEPTH_B",
            "SYMBOL_B",
            "GENEID_B",
            "vep_SV_overlap_name_B",
            "vep_SV_overlap_AF_B",
            "DEL_SYMBOLS",
            "DUP_SYMBOLS",
            "PR",
            "SR",
            "Rescue",
        ],
    ]

    return df


def onco_ts_overlap(s: str, oncogenes_and_ts: set[str]) -> bool:
    """
    Check if any genes in a comma-separated list overlap with oncogenes/TSGs.

    :param s: Comma-separated string of gene symbols
    :param oncogenes_and_ts: Set of oncogene and tumor suppressor gene symbols
    :returns: True if any genes overlap with the oncogene/TSG set
    """

    return not oncogenes_and_ts.isdisjoint(set(s.split(", ")))


def list_all_pairs(a: str, b: str, cosmic_pairs_sorted: set[tuple[str, str]]) -> bool:
    """
    Check if any gene pairs from two lists match known COSMIC fusion pairs.

    Creates all possible pairs from genes in lists a and b, then checks if any match the
    provided set of COSMIC fusion gene pairs.

    :param a: Comma-separated string of gene symbols from breakpoint A
    :param b: Comma-separated string of gene symbols from breakpoint B
    :param cosmic_pairs_sorted: Set of sorted tuples representing known COSMIC fusion
        gene pairs
    :returns: True if any gene pair matches a COSMIC fusion pair
    """

    alist = a.split(", ")
    blist = b.split(", ")

    all_pairs = list(itertools.product(alist, blist))
    all_pairs = set([tuple(sorted(elem)) for elem in all_pairs])

    return not cosmic_pairs_sorted.isdisjoint(all_pairs)


config = {}

with open("./config.toml", "rb") as f:
    config.update(tomllib.load(f))

legacy_workspace = TerraWorkspace(
    workspace_namespace=config["terra"]["legacy_workspace_namespace"],
    workspace_name=config["terra"]["legacy_workspace_name"],
    owners=[],
)

workspace = TerraWorkspace(
    workspace_namespace=config["terra"]["workspace_namespace"],
    workspace_name=config["terra"]["workspace_name"],
    owners=[],
)

src_bucket_name = call_firecloud_api(
    firecloud_api.get_workspace,
    namespace=legacy_workspace.workspace_namespace,
    workspace=legacy_workspace.workspace_name,
    fields=["workspace.bucketName"],
)["workspace"]["bucketName"]

dest_bucket_name = call_firecloud_api(
    firecloud_api.get_workspace,
    namespace=workspace.workspace_namespace,
    workspace=workspace.workspace_name,
    fields=["workspace.bucketName"],
)["workspace"]["bucketName"]

storage_client = storage.Client()
src_bucket = storage_client.bucket(src_bucket_name)
dest_bucket = storage_client.bucket(dest_bucket_name)

samples = pd.read_parquet("~/Desktop/src_samples.parquet")
samples = samples.rename(columns={"sv_selected_somatic_sv": "bedpe"})

samples["sv_selected_somatic_sv_fixed"] = (
    samples["bedpe"]
    .str.replace(".svs.expanded.reannotated.bedpe", ".selected_somatic_sv.parquet")
    .str.replace(src_bucket_name, dest_bucket_name, regex=False)
)

blobs = list_blobs(dest_bucket_name, glob="**/*.selected_somatic_sv.parquet")
samples["done"] = (
    samples["sv_selected_somatic_sv_fixed"].isin(blobs["url"])
    if len(blobs) > 0
    else False
)
samples = samples.sort_values("sample_id")


def fix(r):
    print(r["sample_id"])
    df = pd.read_csv(r["bedpe"], sep="\t", low_memory=False)

    df = filter_svs(
        df,
        cosmic_fusion_gene_pairs_path=Path(
            "/Users/dmccabe/Data/ref/sv/cosmic_fusion_gene_pairs_2025-06-03.csv"
        ),
        onco_tsg_path=Path("/Users/dmccabe/Data/ref/sv/onco_tsg.tsv"),
        min_depth=5,
        sv_gnomad_cutoff=0.001,
        large_sv_size=1000000000,
    )

    df.to_parquet(r["sv_selected_somatic_sv_fixed"])


with ThreadPoolExecutor(max_workers=10) as executor:
    futures = [
        executor.submit(fix, r) for _, r in samples.loc[~samples["done"]].iterrows()
    ]

    for f in futures:
        f.result()  # to raise exceptions, if any

workspace.upload_entities(
    samples.drop(columns=["bedpe", "done"]).rename(
        columns={"sv_selected_somatic_sv_fixed": "sv_selected_somatic_sv"}
    )
)
