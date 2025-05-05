"""
Download the segmental duplications track from UCSC (via Galaxy), then filter and merge
the BED file.
"""

import datetime
import os

import pandas as pd
from bioblend.galaxy import GalaxyInstance
from dotenv import load_dotenv
from natsort import natsort_key

load_dotenv()

today = datetime.datetime.today().strftime("%Y-%m-%d")

gi = GalaxyInstance(url="https://usegalaxy.org", key=os.environ["GALAXY_API_KEY"])
dataset = gi.datasets.get_datasets(name="UCSC Main on Human: genomicSuperDups (genome)")

gi.datasets.download_dataset(
    dataset_id=dataset[0]["id"],
    file_path=f"./data/segdup/genomicSuperDups-{today}.tsv",
    use_default_filename=False,
)

segdup = pd.read_csv(
    f"./data/segdup/genomicSuperDups-{today}.tsv",
    sep="\t",
    low_memory=False,
    usecols=["chrom", "chromStart", "chromEnd", "fracMatch"],
)

segdup = segdup.dropna().astype(
    {
        "chrom": "string",
        "chromStart": "int64",
        "chromEnd": "int64",
        "fracMatch": "float64",
    }
)

segdup = segdup.loc[
    segdup["chrom"].str.match(r"^chr[0-9 X Y]+$") & segdup["fracMatch"].ge(0.98)
].drop(columns="fracMatch")

segdup = segdup.sort_values(by=["chrom", "chromStart", "chromEnd"], key=natsort_key)

segdup["is"] = 1

segdup.to_csv(
    f"./data/segdup/segdup_unmerged-{today}.bed", sep="\t", index=False, header=False
)

"""
bedtools merge -c 4 -o distinct -i data/segdup/segdup_unmerged-2025-05-05.bed \
    > data/segdup/segdup-2025-05-05.bed
bgzip data/segdup/segdup-2025-05-05.bed -k -f
tabix data/segdup/segdup-2025-05-05.bed.gz -0 -s1 -b2 -e3 -f
"""
