"""
Compare the filter annotations in two `FilterMutectCalls` output VCFs for the same
sample to determine whether https://gatk.broadinstitute.org/hc/en-us/community/posts/4404184803227
is fixed, meaning we can remove the `fix_mutect2` task.
"""

import re

import pandas as pd


def read_fp_vcf(path: str) -> pd.DataFrame:
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=[
            "chr",
            "pos",
            "id",
            "ref",
            "alt",
            "qual",
            "filter",
            "info",
            "format",
            "values",
        ],
        usecols=["chr", "pos", "ref", "alt", "filter"],
        comment="#",
    )

    return df


stable = read_fp_vcf("./data/fix_mutect2/stable.vcf")
nightly = read_fp_vcf("./data/fix_mutect2/nightly.vcf")

comp = stable.merge(
    nightly,
    how="inner",
    on=["chr", "pos", "ref", "alt"],
    suffixes=("_stable", "_nightly"),
)

diff = comp.loc[comp["filter_stable"].ne(comp["filter_nightly"])]
filter_diff_counts = (
    diff[["filter_stable", "filter_nightly"]].value_counts().reset_index()
)
filter_diff_counts["prop"] = (filter_diff_counts["count"] / len(comp)).round(3)

removable_filters = [
    "weak_evidence",
    "map_qual",
    "strand_bias",
    "slippage",
    "clustered_events",
    "base_qual",
]

rf_regex = re.compile("|".join(removable_filters))

filter_diff_counts["would_be_removed_stable"] = filter_diff_counts[
    "filter_stable"
].str.contains(rf_regex)

filter_diff_counts["would_be_removed_nightly"] = filter_diff_counts[
    "filter_nightly"
].str.contains(rf_regex)

clustered_either = filter_diff_counts.loc[
    filter_diff_counts["filter_stable"].str.contains("clustered_events")
    | filter_diff_counts["filter_nightly"].str.contains("clustered_events")
]
clustered_both = filter_diff_counts.loc[
    filter_diff_counts["filter_stable"].str.contains("clustered_events")
    & filter_diff_counts["filter_nightly"].str.contains("clustered_events")
]
clustered_only_stable = filter_diff_counts.loc[
    filter_diff_counts["filter_stable"].str.contains("clustered_events")
    & ~filter_diff_counts["filter_nightly"].str.contains("clustered_events")
]
clustered_only_nightly = filter_diff_counts.loc[
    ~filter_diff_counts["filter_stable"].str.contains("clustered_events")
    & filter_diff_counts["filter_nightly"].str.contains("clustered_events")
]

would_be_removed_either = filter_diff_counts.loc[
    filter_diff_counts["would_be_removed_stable"]
    | filter_diff_counts["would_be_removed_nightly"]
]
would_be_removed_both = filter_diff_counts.loc[
    filter_diff_counts["would_be_removed_stable"]
    & filter_diff_counts["would_be_removed_nightly"]
]
would_be_removed_only_stable = filter_diff_counts.loc[
    filter_diff_counts["would_be_removed_stable"]
    & ~filter_diff_counts["would_be_removed_nightly"]
]
would_be_removed_only_nightly = filter_diff_counts.loc[
    ~filter_diff_counts["would_be_removed_stable"]
    & filter_diff_counts["would_be_removed_nightly"]
]

print(
    would_be_removed_only_stable[["filter_stable", "filter_nightly", "count", "prop"]]
)
