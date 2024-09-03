"""
Download the most recent hg38 data set available from (e.g.)
https://cancer.sanger.ac.uk/cosmic/download/cosmic/v99/fusion
"""

import datetime
import re
from collections import OrderedDict

import pandas as pd

VERSION = "99"

fusions = (
    pd.read_csv(
        f"./data/Cosmic_Fusion_v{VERSION}_GRCh38.tsv.gz",
        sep="\t",
        low_memory=False,
        usecols=["FUSION_SYNTAX"],
        dtype="string",
    )["FUSION_SYNTAX"]
    .dropna()
    .drop_duplicates()
)

fusion_dict = OrderedDict()

# adapted from createCosmicFusionGeneTsv.py in https://github.com/broadinstitute/gatk
for x in fusions:
    genes_in_this_fusion = re.findall(r"\(([A-Z0-9-.]+)\)", x)

    for k in genes_in_this_fusion:
        if k not in fusion_dict.keys():
            fusion_dict[k] = dict()

        # Look for the fusion
        if x not in fusion_dict[k].keys():
            fusion_dict[k][x] = 0

        fusion_dict[k][x] = fusion_dict[k][x] + 1


def fmt_fusion(d):
    arr = []
    for k in d.keys():
        s = "%(k)s(%(fgene)s)" % {"k": k, "fgene": str(d[k])}
        arr.append(s)
    return "|".join(arr)


fusions_df = pd.DataFrame(
    {
        "gene": fusion_dict.keys(),
        "fusion_genes": [fmt_fusion(x) for x in fusion_dict.values()],
    }
)

today = datetime.datetime.today().strftime("%Y-%m-%d")

fusions_df.to_csv(f"./data/cosmic_fusions_{today}.csv", index=False)
