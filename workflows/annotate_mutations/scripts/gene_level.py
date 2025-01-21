"""
This script uses the OncoKB API to dopwnload lists of oncogenes and tumor suppressor
genes. The files should then be uploaded to GCS.
"""

import datetime
import os
from io import StringIO

import gffutils
import pandas as pd
import requests
from dotenv import load_dotenv

load_dotenv()

today = datetime.datetime.today().strftime("%Y-%m-%d")

r = requests.get(
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.primary_assembly.annotation.gtf.gz"
)

genecode_gtf_path = "./data/gencode.gtf.gz"

with open(genecode_gtf_path, "wb") as f:
    for chunk in r.iter_content(chunk_size=8192):
        f.write(chunk)

gencode_db_path = "./data/gencode.db"

if not os.path.exists(gencode_db_path):
    gencode_db = gffutils.create_db(
        genecode_gtf_path,
        dbfn=gencode_db_path,
        force=True,
        keep_order=True,
        merge_strategy="merge",
        sort_attribute_values=True,
        disable_infer_genes=True,
        disable_infer_transcripts=True,
    )
else:
    print("Connecting to existing gffutils database...")
    gencode_db = gffutils.FeatureDB(gencode_db_path)

valid_chroms = [*[f"chr{x}" for x in range(1, 23)], "chrX", "chrY"]
genes = []

for feature in gencode_db.features_of_type("gene"):
    chrom = feature.seqid
    start = feature.start
    end = feature.end

    hugo_symbol = feature.attributes.get("gene_name", ["Unknown"])[0]

    if chrom in valid_chroms:
        genes.append([chrom, start, end, hugo_symbol])

gencode_df = pd.DataFrame(genes, columns=["chrom", "start", "end", "hugo_symbol"])

r = requests.get(
    "https://www.oncokb.org/api/v1/utils/cancerGeneList",
    headers={"Authorization": f"Bearer {os.environ['ONCOKB_API_KEY']}"},
)

j = r.json()

onco_tsg = pd.DataFrame(j).astype(
    dtype={
        "hugoSymbol": "string",
        "entrezGeneId": "int64",
        "grch37Isoform": "string",
        "grch37RefSeq": "string",
        "grch38Isoform": "string",
        "grch38RefSeq": "string",
        "oncokbAnnotated": "boolean",
        "occurrenceCount": "int64",
        "mSKImpact": "boolean",
        "mSKHeme": "boolean",
        "foundation": "boolean",
        "foundationHeme": "boolean",
        "vogelstein": "boolean",
        "sangerCGC": "boolean",
        "geneAliases": "object",
        "tsg": "boolean",
        "oncogene": "boolean",
    }
)

onco_tsg = onco_tsg.loc[onco_tsg["oncogene"] | onco_tsg["tsg"]]

onco_tsg = pd.concat(
    [
        onco_tsg.drop(columns="geneAliases"),
        onco_tsg.drop(columns="hugoSymbol")
        .explode("geneAliases")
        .rename(columns={"geneAliases": "hugoSymbol"}),
    ],
    ignore_index=True,
)

onco_tsg_bed = gencode_df.merge(
    onco_tsg,
    how="inner",
    left_on="hugo_symbol",
    right_on="hugoSymbol",
).drop(columns="hugoSymbol")

onco_tsg_bed["chrom"] = onco_tsg_bed["chrom"].str.lstrip("chr")
onco_tsg_bed["chrom"] = (
    onco_tsg_bed["chrom"].replace({"X": "23", "Y": "24"}).astype("int8")
)
onco_tsg_bed = onco_tsg_bed.sort_values(["chrom", "start", "end"])
onco_tsg_bed["chrom"] = (
    onco_tsg_bed["chrom"].astype("string").replace({"23": "X", "24": "Y"})
)
onco_tsg_bed["chrom"] = "chr" + onco_tsg_bed["chrom"]

onco_tsg_bed["chrom_lead"] = onco_tsg_bed["chrom"].shift(-1)
onco_tsg_bed["chrom_lag"] = onco_tsg_bed["chrom"].shift(1)
onco_tsg_bed["start_lead"] = onco_tsg_bed["start"].shift(-1).astype("Int64")
onco_tsg_bed["end_lag"] = onco_tsg_bed["end"].shift(1).astype("Int64")

onco_tsg_bed["overlap"] = False

onco_tsg_bed.loc[
    (
        onco_tsg_bed["chrom"].eq(onco_tsg_bed["chrom_lead"])
        & onco_tsg_bed["start_lead"].between(onco_tsg_bed["start"], onco_tsg_bed["end"])
    )
    | (
        onco_tsg_bed["chrom"].eq(onco_tsg_bed["chrom_lag"])
        & onco_tsg_bed["end_lag"].between(onco_tsg_bed["start"], onco_tsg_bed["end"])
    ),
    "overlap",
] = True

print(onco_tsg_bed.loc[onco_tsg_bed["overlap"]])

onco_tsg_bed["start"] = onco_tsg_bed["start"] - 1
onco_tsg_bed = onco_tsg_bed[["chrom", "start", "end", "oncogene", "tsg"]]
onco_tsg_bed[["oncogene", "tsg"]] = onco_tsg_bed[["oncogene", "tsg"]].astype("int8")
onco_tsg_bed.to_csv(
    f"./data/oncokb/onco_tsg_{today}.bed", header=False, index=False, sep="\t"
)

"""
echo '##INFO=<ID=ONCOGENE,Number=0,Type=Flag,Description="Is oncogene">' \
    > ./data/oncokb/onco_tsg.hdr.vcf
echo '##INFO=<ID=TSG,Number=0,Type=Flag,Description="Is tumor suppressor gene">' \
    >> ./data/oncokb/onco_tsg.hdr.vcf

bgzip data/oncokb/onco_tsg_2024-11-25.bed -k -f
tabix ./data/oncokb/onco_tsg_2024-11-25.bed.gz -0 -s1 -b2 -e3 -f

# e.g.
bcftools annotate ~/Desktop/the.vcf.gz \
    --annotations=./data/oncokb/onco_tsg_2024-11-25.bed.gz \
    --output=./data/oncokb/onco_tsg_annot.vcf \
    --header-lines=./data/oncokb/onco_tsg.hdr.vcf \
    --columns=CHROM,BEG,END,ONCOGENE,TSG \
    --merge-logic="ONCOGENE:first,TSG:first"
"""

r = requests.get(
    "https://www.genenames.org/cgi-bin/download/custom",
    params={
        "col": [
            "gd_app_sym",
            "gd_aliases",
            "gd_app_name",
            "family.name",
        ],
        "status": "Approved",
        "hgnc_dbtag": "on",
        "format": "text",
        "submit": "submit",
    },
)

hgnc = pd.read_csv(StringIO(r.text), sep="\t", dtype="string")

hgnc = hgnc.rename(
    columns={
        "Approved symbol": "hugo_symbol",
        "Alias symbols": "aliases",
        "Approved name": "approved_name",
        "Gene group name": "gene_group_name",
    }
)

hgnc["aliases"] = hgnc["aliases"].str.split(", ", expand=False)

hgnc = pd.concat(
    [
        hgnc.drop(columns="aliases").assign(main_symbol=True),
        hgnc.dropna(subset="aliases")
        .assign(main_symbol=False)
        .drop(columns="hugo_symbol")
        .explode("aliases")
        .rename(columns={"aliases": "hugo_symbol"}),
    ],
    ignore_index=True,
).drop_duplicates()

hgnc_bed = gencode_df.merge(hgnc, how="inner", on="hugo_symbol")

hgnc_bed_main = hgnc_bed.loc[hgnc_bed["main_symbol"]]
hgnc_bed_alias = hgnc_bed.loc[
    ~hgnc_bed["main_symbol"]
    & ~hgnc_bed["hugo_symbol"].isin(hgnc_bed_main["hugo_symbol"])
]

hgnc_bed = pd.concat([hgnc_bed_main, hgnc_bed_alias]).drop(columns="main_symbol")

hgnc_bed["chrom"] = hgnc_bed["chrom"].str.lstrip("chr")
hgnc_bed["chrom"] = hgnc_bed["chrom"].replace({"X": "23", "Y": "24"}).astype("int8")
hgnc_bed = hgnc_bed.sort_values(["chrom", "start", "end"])
hgnc_bed["chrom"] = hgnc_bed["chrom"].astype("string").replace({"23": "X", "24": "Y"})
hgnc_bed["chrom"] = "chr" + hgnc_bed["chrom"]

hgnc_bed["start"] = hgnc_bed["start"] - 1

hgnc_bed["approved_name"] = hgnc_bed["approved_name"].str.replace(",", "%2C")
hgnc_bed["gene_group_name"] = hgnc_bed["gene_group_name"].str.replace(",", "%2C")

hgnc_bed[["chrom", "start", "end", "approved_name", "gene_group_name"]].to_csv(
    f"./data/hgnc_{today}.bed", header=False, index=False, sep="\t"
)

"""
echo '##INFO=<ID=HGNC_NAME,Number=.,Type=String,Description="HGNC approved name">' \
    > ./data/hgnc.hdr.vcf
echo '##INFO=<ID=HGNC_GROUP,Number=.,Type=String,Description="HGNC gene group name">' \
    >> ./data/hgnc.hdr.vcf

bgzip data/hgnc_2024-11-25.bed -k -f
tabix ./data/hgnc_2024-11-25.bed.gz -0 -s1 -b2 -e3 -f

# e.g.
bcftools annotate ~/Desktop/the.vcf.gz \
    --annotations=./data/hgnc_2024-11-25.bed.gz \
    --output=./data/hgnc_annot.vcf \
    --header-lines=./data/hgnc.hdr.vcf \
    --columns=CHROM,BEG,END,HGNC_NAME,HGNC_GROUP \
    --merge-logic="HGNC_NAME:unique,HGNC_GROUP:unique"
"""
