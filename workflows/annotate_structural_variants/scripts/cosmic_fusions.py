"""
Download the most recent hg38 data set available from (e.g.)
https://cancer.sanger.ac.uk/cosmic/download/cosmic/v99/fusion
"""

import datetime
import re
from collections import OrderedDict

import pandas as pd

VERSION = "102"

fusions = (
    pd.read_csv(
        f"./data/Cosmic_Fusion_v{VERSION}_GRCh38.tsv.gz",
        sep="\t",
        low_memory=False,
        usecols=["FIVE_PRIME_GENE_SYMBOL", "THREE_PRIME_GENE_SYMBOL"],
        dtype="string",
    )
    .dropna()
    .drop_duplicates()
)

fusions.columns = fusions.columns.str.lower()

today = datetime.datetime.today().strftime("%Y-%m-%d")

fusions.to_csv(f"./data/cosmic_fusion_gene_pairs_{today}.csv", index=False)
