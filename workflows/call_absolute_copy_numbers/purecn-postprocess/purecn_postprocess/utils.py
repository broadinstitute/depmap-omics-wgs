import json
from pathlib import Path
from typing import Literal, Type

import pandas as pd
from click import echo

from purecn_postprocess.types import Loh, PanderaBaseSchema, Solution, TypedDataFrame


def collect_outputs(solution_path: Path, loh_path: Path, out_path: Path) -> None:
    """
    Collect outputs from PureCN solution and gene call data, call chromosomal
    instability (CIN) and whole genome doubling (WGD), then output the results as JSON.

    :param solution_path: path to the PureCN solution CSV file
    :param loh_path: path to the `_loh.csv` (gene calls) CSV file
    :param out_path: path to the output JSON file
    """

    # read PureCN output CSVs
    sol = read_solution(solution_path)
    loh = read_loh(loh_path)

    # call chromosomal instability (CIN)
    cin = calculate_cin(loh, allele_specific=False, reference_state="normal")
    cin_ploidy_robust = calculate_cin(
        loh, allele_specific=False, reference_state="dominant"
    )
    cin_allele_specific = calculate_cin(
        loh, allele_specific=True, reference_state="normal"
    )
    cin_allele_specific_ploidy_robust = calculate_cin(
        loh, allele_specific=True, reference_state="dominant"
    )

    # call whole genome doubling (WGD)
    is_wgd = call_wgd(loh, sol["ploidy"].squeeze())

    # output everything as JSON to read back in inside the workflow
    res = {
        **sol.to_dict(orient="records")[0],
        "cin": cin,
        "cin_ploidy_robust": cin_ploidy_robust,
        "cin_allele_specific": cin_allele_specific,
        "cin_allele_specific_ploidy_robust": cin_allele_specific_ploidy_robust,
        "is_wgd": is_wgd,
    }

    # write JSON to file
    with open(out_path, "w") as f:
        f.write(json.dumps(res) + "\n")


def read_solution(solution_path: Path) -> TypedDataFrame[Solution]:
    """
    Read the PureCN solution CSV file and return a data frame with the relevant columns.

    :param solution_path: path to the PureCN solution CSV file
    :return: data frame containing the relevant columns from the PureCN solution
    """

    echo(f"Reading PureCN solution from {solution_path}")

    purecn_df = pd.read_csv(
        solution_path,
        dtype={
            "Sampleid": "string",
            "Purity": "float64",
            "Ploidy": "float64",
            "Sex": pd.CategoricalDtype(categories=["M", "F"]),
            "Contamination": "float64",
            "Flagged": "boolean",
            "Failed": "boolean",
            "Curated": "boolean",
            "Comment": "string",
        },
        usecols=["Purity", "Ploidy", "Contamination", "Flagged", "Curated", "Comment"],
    )

    assert len(purecn_df) == 1
    purecn_df.columns = purecn_df.columns.str.lower()
    return type_data_frame(purecn_df, Solution)


def read_loh(loh_path: Path) -> TypedDataFrame[Loh]:
    """
    Read the PureCN gene calls CSV file and return a data frame with the relevant
    columns.

    :param loh_path: path to the `_loh.csv` (gene calls) CSV file
    :return: data frame containing the relevant columns from the PureCN gene call data
    """

    echo(f"Reading gene calls from {loh_path}")

    loh = pd.read_csv(
        loh_path,
        dtype={
            "Sampleid": "string",
            "chr": "string",
            "start": "int64",
            "end": "int64",
            "arm": pd.CategoricalDtype(categories=["p", "q"]),
            "C": "float64",
            "M": "Int64",
            "type": "string",
            "seg.mean": "float64",
            "num.mark": "int64",
            "num.snps": "int64",
            "M.flagged": "boolean",
            "maf.expected": "float64",
            "maf.observed": "float64",
        },
        usecols=["start", "end", "C", "M", "type"],
    )

    loh.columns = loh.columns.str.lower()

    loh["size"] = loh["end"] - loh["start"] + 1
    loh = loh.dropna(subset=["size"])

    # PureCN outputs non-integer segment copy numbers when they exceed `max.copy.number`
    loh["c"] = loh["c"].round(0).astype("int64")

    return type_data_frame(loh, Loh)


def calculate_cin(
    loh: TypedDataFrame[Loh],
    allele_specific: bool,
    reference_state: Literal["dominant", "normal"],
) -> float:
    """
    Calculate chromosomal instability (CIN).

    :param loh: data frame containing the gene call data
    :param allele_specific: whether to consider allele-specific copy number variations
    :param reference_state: the reference state to compare against (domainant or normal)
    :return: the proportion of the genome matching the reference state
    """

    echo(
        "Calling chromosomal instability given"
        f" allele-specific={allele_specific} and reference_state='{reference_state}'"
    )

    loh = loh.copy()

    if allele_specific:
        # remove rows with NaN in 'M' (minor integer copy number) column
        loh = loh.dropna(subset=["m"])

        # write as, e.g., "2/1"
        loh["state"] = loh["c"].astype("str") + "/" + loh["m"].astype("str")

    else:
        loh["state"] = loh["c"].astype("str")

    # get the index name ("state") of dominant copy number state
    reference_state_cn = loh.groupby("state")["size"].sum().idxmax()

    if reference_state == "normal":
        reference_state_cn = "2/1" if allele_specific else "2"

    # calculate proportion of genome not matching the reference state
    loh["is_reference"] = loh["state"].eq(reference_state_cn)
    return loh.loc[~loh["is_reference"], "size"].sum() / loh["size"].sum()


def call_wgd(loh: TypedDataFrame[Loh], ploidy: float) -> bool:
    """
    Call whole genome doubling (WGD). Citation:
    https://www.sciencedirect.com/science/article/pii/S0092867421002944

    :param loh: data frame containing the gene call data
    :param ploidy: the ploidy of the sample from the solution data
    :return: whether the sample has undergone WGD
    """

    echo(f"Calling whole genome doubling given ploidy={ploidy}")

    loh = loh.copy()

    # use LOH regions to calculate proportion of genome having undergone WGD
    loh_len = loh.loc[
        loh["type"].isin({"LOH", "COPY-NEUTRAL LOH", "WHOLE ARM COPY-NEUTRAL LOH"}),
        "size",
    ].sum()

    loh_frac = loh_len / loh["size"].sum()
    stat = -2 * loh_frac + 3
    return bool(stat < ploidy)


def type_data_frame(
    df: pd.DataFrame,
    pandera_schema: Type[PanderaBaseSchema],
    remove_unknown_cols: bool = False,
) -> TypedDataFrame[PanderaBaseSchema]:
    """
    Coerce a data frame into one specified by a Pandera schema and optionally remove
    unknown columns.

    :param df: a data frame
    :param pandera_schema: a Pandera schema
    :param remove_unknown_cols: remove columns not specified in the schema
    :return: a data frame validated with the provided Pandera schema
    """

    if len(df) == 0:
        # make an empty data frame that conforms to the Pandera schema
        s = pandera_schema.to_schema()

        # `example` doesn't know how to instantiate dicts, so do that manually
        dict_cols = []

        for c in s.columns:
            if s.columns[c].dtype.type is dict:
                dict_cols.append(c)
                s = s.remove_columns([c])

        df = pd.DataFrame(s.example(size=0))

        if len(dict_cols) > 0:
            for c in dict_cols:
                df[c] = {}

    elif remove_unknown_cols:
        df_cols = pandera_schema.to_schema().columns.keys()
        df = df.loc[:, df_cols]

    return TypedDataFrame[pandera_schema](df)
