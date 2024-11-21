from io import StringIO

import pandas as pd
import requests

r = requests.get(
    "https://www.genenames.org/cgi-bin/download/custom",
    params={
        "col": ["gd_hgnc_id", "gd_app_name", "family.name"],
        "status": "Approved",
        "hgnc_dbtag": "on",
        "order_by": "gd_app_sym_sort",
        "format": "text",
        "submit": "submit",
    },
)

df = pd.read_csv(StringIO(r.text), sep="\t", dtype="string")

df = df.rename(
    columns={
        "HGNC ID": "hgnc_id",
        "Approved name": "hgnc_name",
        "Gene group name": "hgnc_family",
    }
)

df.to_csv("./data/hgnc.tsv", sep="\t", index=False)
