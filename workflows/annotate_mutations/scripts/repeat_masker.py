"""
Download the RepeatMAsker track from UCSC (via Galaxy), then filter and merge the BED
file.
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
dataset = gi.datasets.get_datasets(name="UCSC Main on Human: rmsk (genome)")

gi.datasets.download_dataset(
    dataset_id=dataset[0]["id"],
    file_path=f"./data/repeat_masker/repeatMasker-{today}.tsv",
    use_default_filename=False,
)

repeat_masker = pd.read_csv(
    f"./data/repeat_masker/repeatMasker-{today}.tsv",
    sep="\t",
    low_memory=False,
    usecols=["milliDiv", "milliDel", "milliIns", "genoName", "genoStart", "genoEnd"],
)

repeat_masker = repeat_masker.dropna().astype(
    {
        "milliDiv": "float64",
        "milliDel": "float64",
        "milliIns": "float64",
        "genoName": "string",
        "genoStart": "int64",
        "genoEnd": "int64",
    }
)

repeat_masker["max_repeats"] = repeat_masker[["milliDiv", "milliDel", "milliIns"]].max(
    axis=1
)


repeat_masker = repeat_masker.loc[
    repeat_masker["genoName"].str.match(r"^chr[0-9 X Y]+$")
    & repeat_masker["max_repeats"].le(10)
].drop(columns=["max_repeats", "milliDiv", "milliDel", "milliIns"])

repeat_masker = repeat_masker.sort_values(
    by=["genoName", "genoStart", "genoEnd"], key=natsort_key
)

repeat_masker["is"] = 1

repeat_masker.to_csv(
    f"./data/repeat_masker/repeat_masker_unmerged-{today}.bed",
    sep="\t",
    index=False,
    header=False,
)

"""
bedtools merge -c 4 -o distinct -i data/repeat_masker/repeat_masker_unmerged-2025-05-05.bed \
    > data/repeat_masker/repeat_masker-2025-05-05.bed
bgzip data/repeat_masker/repeat_masker-2025-05-05.bed -k -f
tabix data/repeat_masker/repeat_masker-2025-05-05.bed.gz -0 -s1 -b2 -e3 -f
"""
