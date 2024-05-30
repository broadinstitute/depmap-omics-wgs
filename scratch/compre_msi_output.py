import os
from functools import reduce

import numpy as np
import pandas as pd
from firecloud import api as firecloud_api
from google.cloud import storage

from omics_wgs_pipeline.utils import expand_dict_columns

pd.set_option("display.max_columns", 30)
pd.set_option("display.max_colwidth", 50)
pd.set_option("display.max_info_columns", 30)
pd.set_option("display.max_info_rows", 20)
pd.set_option("display.max_rows", 20)
pd.set_option("display.max_seq_items", None)
pd.set_option("display.width", 200)
pd.set_option("expand_frame_repr", True)
pd.set_option("mode.chained_assignment", "warn")

# get the absolute path to the root directory of the repo
root_dir = os.path.abspath(os.path.join(os.path.dirname(__name__)))


def get_terra_outputs(namespace, workspace):
    wgs_cn_samples = pd.DataFrame(
        firecloud_api.get_entities(
            namespace=namespace, workspace=workspace, etype="sample"
        ).json()
    )
    wgs_cn_samples = expand_dict_columns(wgs_cn_samples, name_columns_with_parent=False)
    wgs_cn_samples = wgs_cn_samples.convert_dtypes()
    wgs_cn_samples = wgs_cn_samples.rename(columns={"name": "sample_id"})

    return (
        wgs_cn_samples[
            [
                "sample_id",
                "msisensor2_score",
                "msisensor2_output",
                "msisensor2_output_dis",
                "msisensor2_output_somatic",
            ]
        ]
        .dropna()
        .astype(
            {
                "sample_id": "string",
                "msisensor2_score": "float64",
                "msisensor2_output": "string",
                "msisensor2_output_dis": "string",
                "msisensor2_output_somatic": "string",
            }
        )
    )


def save_msi_outputs(storage_client, df, old_or_new):
    local_dir = os.path.join(root_dir, "data", old_or_new)

    for i in range(len(df)):
        for c in [
            "msisensor2_output",
            "msisensor2_output_dis",
            "msisensor2_output_somatic",
        ]:
            sample_dir = os.path.join(local_dir, df["sample_id"].iloc[i])
            os.makedirs(sample_dir, exist_ok=True)
            local_file = os.path.join(
                sample_dir, ".".join([df["sample_id"].iloc[i], c])
            )

            blob = storage.Blob.from_string(df[c].iloc[i], client=storage_client)
            blob.download_to_filename(local_file)

            df.iloc[i, df.columns.get_loc(c)] = local_file

    return df


def compare_scores(old_samples, new_samples):
    score_comp = old_samples[["sample_id", "msisensor2_score"]].merge(
        new_samples[["sample_id", "msisensor2_score"]],
        on="sample_id",
        how="inner",
        suffixes=("_old", "_new"),
    )

    score_med_abs_dev = np.abs(
        score_comp["msisensor2_score_new"] - score_comp["msisensor2_score_old"]
    ).median()
    score_max_abs_dev = np.abs(
        score_comp["msisensor2_score_new"] - score_comp["msisensor2_score_old"]
    ).max()

    print(f"Median absolute deviation: {round(score_med_abs_dev, 2)}")
    print(f"Max absolute deviation: {round(score_max_abs_dev, 2)}")


def get_weighted_avg(dis):
    dis = pd.Series(dis.replace("T:", "").strip().split(" ")).astype(int)
    return np.ma.average(dis.index + 1, weights=dis.values)


def get_n_repeats_from_msisensor2_output_dis(path):
    df = pd.read_csv(path, sep="\t", header=None)

    # split site rows from distribution rows (they alternate)
    df = pd.DataFrame(
        {
            "site": df.loc[df.index % 2 == 0, 0].values,
            "dis": df.loc[df.index % 2 == 1, 0].values,
        }
    )

    # split the site into its component info
    site_cols = ["chr", "loc", "left_flank", "repeat", "right_flank"]

    df = pd.concat(
        [
            pd.DataFrame(df["site"].str.split(" ").tolist(), columns=site_cols),
            df["dis"],
        ],
        axis=1,
    )

    # calc weighted average number of repeats
    df["wavg"] = df["dis"].apply(get_weighted_avg)
    return df


def combine_dis_wavg(df):
    repeats = [
        get_n_repeats_from_msisensor2_output_dis(x["msisensor2_output_dis"])
        .drop(columns="dis")
        .rename(columns={"wavg": x["sample_id"]})
        .set_index(["chr", "loc", "left_flank", "repeat", "right_flank"])
        for _, x in df.iterrows()
    ]

    # reduce a list of data frames into a single dataframe using `join`
    return reduce(lambda x, y: x.join(y, how="outer"), repeats)


def compare_repeats(old_samples, new_samples):
    old_repeats = combine_dis_wavg(old_samples)
    new_repeats = combine_dis_wavg(new_samples)

    repeats_comp = old_repeats.join(
        new_repeats, how="outer", lsuffix="_old", rsuffix="_new"
    )

    # make sure the outer join was effectively an inner one (i.e. indexes match
    assert len(repeats_comp) == len(old_repeats) == len(new_repeats)

    id_vars = ["chr", "loc", "left_flank", "repeat", "right_flank"]

    old_repeats_long = old_repeats.reset_index().melt(
        id_vars=id_vars,
        var_name="sample_id",
        value_name="wavg",
    )

    new_repeats_long = new_repeats.reset_index().melt(
        id_vars=id_vars,
        var_name="sample_id",
        value_name="wavg",
    )
    long_comp = old_repeats_long.merge(
        new_repeats_long,
        on=[*id_vars, "sample_id"],
        how="inner",
        suffixes=("_old", "_new"),
    )

    long_comp["abs_dev"] = np.abs(long_comp["wavg_old"] - long_comp["wavg_new"])
    overall_med_abs_dev = long_comp["abs_dev"].median()
    print(f"Overall median absolute deviation: {round(overall_med_abs_dev, 2)}")

    long_comp["abs_dev_rank"] = long_comp["abs_dev"].rank(ascending=False, method="min")


old_samples = get_terra_outputs(
    namespace="broad-firecloud-ccle", workspace="DepMap_WGS_CN"
)
new_samples = get_terra_outputs(
    namespace="broad-firecloud-ccle", workspace="omics_wgs_pipeline"
)

sample_ids = set(old_samples["sample_id"]).intersection(set(new_samples["sample_id"]))
old_samples = old_samples.loc[old_samples["sample_id"].isin(sample_ids)]
new_samples = new_samples.loc[new_samples["sample_id"].isin(sample_ids)]

storage_client = storage.Client(project="depmap-omics")
old_samples = save_msi_outputs(storage_client, old_samples, old_or_new="old")
new_samples = save_msi_outputs(storage_client, new_samples, old_or_new="new")

compare_scores(old_samples, new_samples)
