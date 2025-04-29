"""
This script downloads the latest data form CivicDB and transforms the data to be
uploaded to GCs.
"""

import datetime
from urllib.parse import quote

import pandas as pd
import requests
from natsort import natsort_key
from pyliftover import LiftOver

liftover = LiftOver("hg19", "hg38")


def liftover_one(chrom, pos):
    lifted = liftover.convert_coordinate(chrom, pos)

    if len(lifted) != 1:
        raise NotImplementedError

    return ",".join([lifted[0][0], str(lifted[0][1])])


today = datetime.datetime.today().strftime("%Y-%m-%d")


with open(f"./data/civic/variant_summaries_{today}.tsv", "w") as f:
    r = requests.get(
        "https://civicdb.org/downloads/nightly/nightly-VariantSummaries.tsv"
    )
    f.write(r.text)

with open(f"./data/civic/molecular_profile_summaries_{today}.tsv", "w") as f:
    r = requests.get(
        "https://civicdb.org/downloads/nightly/nightly-MolecularProfileSummaries.tsv"
    )
    f.write(r.text)

vs = pd.read_csv(
    f"./data/civic/variant_summaries_{today}.tsv",
    sep="\t",
    dtype={
        "variant_id": "string",
        "variant_civic_url": "string",
        "gene": "string",
        "entrez_id": "Int64",
        "variant": "string",
        "variant_groups": "string",
        "chromosome": "string",
        "start": "Int64",
        "stop": "Int64",
        "reference_bases": "string",
        "variant_bases": "string",
        "representative_transcript": "string",
        "ensembl_version": "string",
        "reference_build": "string",
        "chromosome2": "string",
        "start2": "Int64",
        "stop2": "Int64",
        "representative_transcript2": "string",
        "variant_types": "string",
        "hgvs_descriptions": "string",
        "last_review_date": "string",
        "allele_registry_id": "string",
        "clinvar_ids": "string",
        "variant_aliases": "string",
        "is_flagged": "boolean",
        "single_variant_molecular_profile_id": "int64",
    },
)

mps = pd.read_csv(
    f"./data/civic/molecular_profile_summaries_{today}.tsv",
    sep="\t",
    dtype={
        "name": "string",
        "molecular_profile_id": "int64",
        "summary": "string",
        "variant_ids": "string",
        "variants_civic_url": "string",
        "evidence_score": "float64",
        "evidence_item_ids": "string",
        "evidence_items_civic_url": "string",
        "assertion_ids": "string",
        "assertions_civic_url": "string",
        "aliases": "string",
        "last_review_date": "string",
        "is_flagged": "boolean",
    },
)

mps["variant_ids"] = mps["variant_ids"].str.split(r",\s*")
mps = mps.explode(column="variant_ids", ignore_index=True)

df = vs.merge(mps, how="inner", left_on="variant_id", right_on="variant_ids").rename(
    columns={
        "chromosome": "chrom",
        "start": "pos",
        "reference_bases": "ref",
        "variant_bases": "alt",
        "summary": "description",
    }
)

df["chrom"] = "chr" + df["chrom"]

df = df.dropna(subset=["variant_id", "chrom", "pos", "ref", "alt"])
df = df.loc[~df[["evidence_score", "description"]].isna().all(axis=1)]

df = df.loc[
    :,
    [
        "variant_id",
        "chrom",
        "pos",
        "ref",
        "alt",
        "evidence_score",
        "description",
        "reference_build",
    ],
]

df.loc[df["reference_build"].eq("GRCh37"), "chrom_pos_hg38"] = df.loc[
    df["reference_build"].eq("GRCh37")
].apply(lambda x: liftover_one(x["chrom"], x["pos"] - 1), axis=1)

df.loc[df["reference_build"].eq("GRCh38"), "chrom_pos_hg38"] = (
    df.loc[df["reference_build"].eq("GRCh38"), "chrom"]
    + ","
    + df.loc[df["reference_build"].eq("GRCh38"), "pos"].astype("string")
)

print(df.loc[df["chrom_pos_hg38"].isna()])
df = df.dropna(subset="chrom_pos_hg38")

df[["chrom_hg38", "pos_hg38"]] = df["chrom_pos_hg38"].str.split(",", expand=True)
df["chrom_hg38"] = df["chrom_hg38"].astype("string")
df["pos_hg38"] = df["pos_hg38"].astype("Int64")
df["pos_hg38"] += 1

df["chrom"] = df["chrom_hg38"]
df["pos"] = df["pos_hg38"]
df = df.drop(columns=["chrom_pos_hg38", "chrom_hg38", "pos_hg38"])

df = df.sort_values(by=["chrom", "pos"], key=natsort_key)

df.loc[~df["description"].isna(), "description"] = df.loc[
    ~df["description"].isna(), "description"
].apply(quote)

df[
    ["chrom", "pos", "ref", "alt", "variant_id", "evidence_score", "description"]
].to_csv(f"./data/civic/civic_{today}.tsv", sep="\t", index=False, header=False)

"""
echo '##INFO=<ID=CIVIC_ID,Number=1,Type=Integer,Description="CIVIC variant ID">' \
    > ./data/civic/civic.hdr.vcf
echo '##INFO=<ID=CIVIC_SCORE,Number=1,Type=String,Description="CIVIC evidence score">' \
    >> ./data/civic/civic.hdr.vcf
echo '##INFO=<ID=CIVIC_DESC,Number=1,Type=String,Description="CIVIC variant description">' \
    >> ./data/civic/civic.hdr.vcf
bgzip ./data/civic/civic_2025-04-29.tsv -k -f
tabix ./data/civic/civic_2025-04-29.tsv.gz -s1 -b2 -e2 -f

# e.g.
bcftools annotate input.vcf.gz \
    --annotations=./data/civic/civic_2025-04-29.tsv.gz \
    --output=./data/output.vcf.gz \
    --header-lines=./data/hess.hdr.vcf \
    --columns=CHROM,POS,REF,ALT,CIVIC_ID,CIVIC_SCORE,CIVIC_DESC
"""
