import datetime

import pandas as pd
from _csv import QUOTE_NONE

df = pd.read_table(
    "https://storage.googleapis.com/broad-public-datasets/funcotator/funcotator_dataSources.v1.7.20200521s/dna_repair_genes/hg19/dnaRepairGenes.20180524T145835.csv",
    sep="|",
    quoting=QUOTE_NONE,
    dtype="string",
    on_bad_lines=lambda x: x[:4],
    engine="python",
)

df = df.rename(
    columns={
        "Gene Name": "gene_name",
        "Activity linked to OMIM": "activity",
        "Chromosome location linked to NCBI MapView": "chrom_pos",
        "Accession number linked to NCBI Entrez": "accession_number",
    }
)

today = datetime.datetime.today().strftime("%Y-%m-%d")

df.to_csv(f"./data/dna_repair_genes_{today}.csv", index=False)
