"""
This script uses the OncoKB API to dopwnload lists of oncogens and tumor suppressor
genes. The files should then be uploaded to GCS.
"""

import datetime
import os

import pandas as pd
import requests
from dotenv import load_dotenv

load_dotenv()

r = requests.get(
    "https://www.oncokb.org/api/v1/utils/cancerGeneList",
    headers={"Authorization": f"Bearer {os.environ['ONCOKB_API_KEY']}"},
)

j = r.json()

df = (
    pd.DataFrame(j)
    .astype(
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
            "geneAliases": "string",
            "tsg": "boolean",
            "oncogene": "boolean",
        }
    )
    .sort_values("hugoSymbol")
)

oncogenes = df.loc[df["oncogene"], "hugoSymbol"].drop_duplicates().to_list()
tsg = df.loc[df["tsg"], "hugoSymbol"].drop_duplicates().to_list()

today = datetime.datetime.today().strftime("%Y-%m-%d")

with open(f"./data/oncokb/oncogenes_{today}.txt", "w") as f:
    f.write("\n".join(oncogenes))

with open(f"./data/oncokb/tsg_{today}.txt", "w") as f:
    f.write("\n".join(tsg))
