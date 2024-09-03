"""
Download the most recent hg38 data set available from (e.g.)
https://cancer.sanger.ac.uk/cosmic/download/cosmic/v99/cancergenecensus
"""

import datetime

import pandas as pd

VERSION = "99"

df = (
    pd.read_csv(
        f"./data/Cosmic_CancerGeneCensus_v{VERSION}_GRCh38.tsv.gz",
        sep="\t",
        low_memory=False,
        usecols=["GENE_SYMBOL", "TRANSLOCATION_PARTNER"],
        dtype="string",
        na_values="?",
    )[["GENE_SYMBOL", "TRANSLOCATION_PARTNER"]]
    .dropna()
    .drop_duplicates()
    .rename(
        columns={
            "GENE_SYMBOL": "gene",
            "TRANSLOCATION_PARTNER": "translocation_partners",
        }
    )
)

today = datetime.datetime.today().strftime("%Y-%m-%d")

df.to_csv(f"./data/cosmic_translocation_partners_{today}.csv", index=False)
