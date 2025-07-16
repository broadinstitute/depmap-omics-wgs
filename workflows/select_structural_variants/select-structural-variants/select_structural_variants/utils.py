import itertools
from pathlib import Path
from typing import Any, Dict, Iterable, List, Set, Tuple, Union

import pandas as pd


def do_select_structural_variants(
    input_bedpe_path: Path,
    gene_annotation_path: Path,
    del_annotation_path: Path,
    dup_annotation_path: Path,
    cosmic_fusion_gene_pairs_path: Path,
    onco_tsg_path: Path,
    out_path: Path,
) -> None:
    """
    Select and filter structural variants from a BEDPE file with gene annotations.

    :param input_bedpe_path: Path to input BEDPE file containing structural variants
    :param gene_annotation_path: Path to gene annotation file for breakpoints
    :param del_annotation_path: Path to gene annotation file for deletions
    :param dup_annotation_path: Path to gene annotation file for duplications
    :param cosmic_fusion_gene_pairs_path: Path to COSMIC fusion gene pairs file
    :param onco_tsg_path: Path to oncogene and tumor suppressor gene file
    :param out_path: Path where filtered results will be saved as parquet
    """

    df = bedpe_to_df(input_bedpe_path)

    df = reannotate_genes(
        df, gene_annotation_path, del_annotation_path, dup_annotation_path
    )

    df = filter_svs(df, cosmic_fusion_gene_pairs_path, onco_tsg_path)
    df.to_parquet(out_path, index=False)


def bedpe_to_df(input_bedpe_path: Path) -> pd.DataFrame:
    """
    Transform a BEDPE file into a DataFrame. Parses the BEDPE format including INFO
    fields and sample-specific FORMAT fields.

    :param input_bedpe_path: Filepath to the BEDPE file
    :returns: DataFrame representation of the BEDPE file
    """

    uniqueargs = ["IMPRECISE"]

    vep_csq_desc = "Consequence annotations from Ensembl VEP."

    with open(input_bedpe_path, "r") as f:
        description, colnames, nrows_toskip = read_comments(f, vep_csq_desc)

    colnames = [i for i in colnames]

    df = pd.read_csv(
        input_bedpe_path,
        sep="\t",
        index_col=False,
        header=None,
        names=colnames,
        skiprows=nrows_toskip,
    )

    vep_fields_a = [k + "_A" for k, v in description.items() if vep_csq_desc in v]
    vep_fields_b = [k + "_B" for k, v in description.items() if vep_csq_desc in v]
    fields_a = {k + "_A": [] for k, _ in description.items()}
    fields_b = {k + "_B": [] for k, _ in description.items()}
    fields_combined = {}

    for suffix, side, fields, vep_fields in [
        ("_A", "INFO_A", fields_a, vep_fields_a),
        ("_B", "INFO_B", fields_b, vep_fields_b),
    ]:
        for j, info in enumerate(df[side].str.split(";").values.tolist()):
            res = {}

            for annot in info:
                if annot == ".":
                    pass

                elif annot in uniqueargs:
                    res.update({annot: True})

                elif "=" in annot:
                    # taking care of the funcotator special fields
                    if "CSQ=" in annot:
                        annot = annot.replace("CSQ=", "")
                        res.update({name: [] for name in vep_fields})

                        for site in annot.split(","):
                            for i, sub_annot in enumerate(site.split("|")):
                                res[vep_fields[i]].append(sub_annot)
                        for k in vep_fields:
                            if "".join(res[k]) != "":
                                res[k] = ",".join(res[k])
                            else:
                                res[k] = ""

                    elif "END" not in annot:
                        k, annot = annot.split("=")
                        res.update({k + suffix: annot})
                else:
                    raise ValueError("unknown argument: " + annot)

            for k in list(fields.keys()):
                fields[k].append(res.get(k, None))

        fields_combined.update(fields)

    df = pd.concat(
        [
            df.drop(columns=["INFO_A", "INFO_B"]),
            pd.DataFrame(data=fields_combined, index=df.index),
        ],
        axis=1,
    )

    samples = [i for i in colnames[21:]]
    sorting = df["FORMAT"][0].split(":")

    for sample in samples:
        res = df[sample].str.split(":").values.tolist()
        maxcols = max([len(v) for v in res])

        if maxcols - len(sorting) > 0:
            for i in range(maxcols - len(sorting)):
                sorting.append(sorting[-1] + "_" + str(i + 1))

        if len(samples) > 1:
            sorting = [sample + "_" + v for v in sorting]

        df = pd.concat(
            [
                df.drop(columns=sample),
                pd.DataFrame(data=res, columns=sorting, index=df.index),
            ],
            axis=1,
        )

    return df


def read_comments(
    f: Iterable[Union[str, bytes]], vep_csq_desc: str
) -> Tuple[Dict[str, str], List[str], int]:
    """
    Read and parse header comments from a BEDPE file.

    Extracts metadata including INFO field descriptions and column names from the file
    header comments.

    :param f: File handle or iterable of lines from the BEDPE file
    :param vep_csq_desc: Description string for VEP CSQ annotations
    :returns: Tuple of (description dict, column names list, rows to skip)
    """

    description = {}
    colnames = []
    rows = 0

    for l in f:
        l = l.decode("utf-8") if isinstance(l, bytes) else l

        if l.startswith("##"):
            rows += 1

            if "INFO" in l[:20]:
                res = l.split("ID=")[1].split(",")[0]

                if res == "CSQ":
                    # print("parsing VEP CSQ")
                    for val in l.split("Description=")[1][:-3].split("|"):
                        val = "vep_" + val.split("Format: ")[-1]
                        description.update({val: vep_csq_desc})
                elif res != "END":
                    desc = l.split("Description=")[1][:-2]
                    description.update({res: desc})

        elif l.startswith("#"):
            colnames = l[1:-1].split("\t")
            rows += 1

        else:
            break

    return description, colnames, rows


def reannotate_genes(
    df: pd.DataFrame,
    gene_annotation_path: Path,
    del_annotation_path: Path,
    dup_annotation_path: Path,
) -> pd.DataFrame:
    """
    Re-annotate genes for structural variants since VEP annotations are unreliable.

    Merges gene annotations for breakpoints, deletions, and duplications with the
    structural variant DataFrame. Splits multi-gene annotations into separate
    symbol and gene ID columns.

    :param df: SVs in BEDPE format
    :param gene_annotation_path: Path to gene re-annotation file for breakpoints
    :param del_annotation_path: Path to gene re-annotation file for DELs only,
        considering the entire deleted intervals
    :param dup_annotation_path: Path to gene re-annotation file for DUPs only,
        considering the entire duplicated intervals
    :returns: SVs in BEDPE format with corrected gene annotations
    """

    # read in annotation files for all breakpoints, DELs, and DUPs separately
    gene_annotation = pd.read_csv(
        gene_annotation_path,
        sep="\t",
        names=[
            "CHROM_A",
            "START_A",
            "END_A",
            "NAME_A",
            "GENE_A",
            "CHROM_B",
            "START_B",
            "END_B",
            "NAME_B",
            "GENE_B",
        ],
    ).astype(
        {
            "CHROM_A": "string",
            "START_A": "int64",
            "END_A": "int64",
            "NAME_A": "string",
            "GENE_A": "string",
            "CHROM_B": "string",
            "START_B": "int64",
            "END_B": "int64",
            "NAME_B": "string",
            "GENE_B": "string",
        }
    )

    del_annotation = pd.read_csv(
        del_annotation_path,
        sep="\t",
        names=["CHROM_A", "START_A", "END_B", "NAME_A", "NAME_B", "DELGENES"],
    ).astype(
        {
            "CHROM_A": "string",
            "START_A": "int64",
            "END_B": "int64",
            "NAME_A": "string",
            "NAME_B": "string",
            "DELGENES": "string",
        }
    )

    dup_annotation = pd.read_csv(
        dup_annotation_path,
        sep="\t",
        names=["CHROM_A", "START_A", "END_B", "NAME_A", "NAME_B", "DUPGENES"],
    ).astype(
        {
            "CHROM_A": "string",
            "START_A": "int64",
            "END_B": "int64",
            "NAME_A": "string",
            "NAME_B": "string",
            "DUPGENES": "string",
        }
    )

    # merge annotations to the SV table one by one
    merged = pd.merge(
        df,
        gene_annotation[
            [
                "CHROM_A",
                "START_A",
                "END_A",
                "NAME_A",
                "GENE_A",
                "CHROM_B",
                "START_B",
                "END_B",
                "GENE_B",
            ]
        ],
        on=["CHROM_A", "START_A", "END_A", "NAME_A", "CHROM_B", "START_B", "END_B"],
        how="outer",
        indicator=True,
    )

    assert merged["_merge"].eq("both").all()

    merged = pd.merge(
        merged.drop(columns="_merge"),
        del_annotation[["CHROM_A", "START_A", "END_B", "NAME_A", "DELGENES"]],
        on=["CHROM_A", "START_A", "END_B", "NAME_A"],
        how="left",
    )

    merged = pd.merge(
        merged,
        dup_annotation[["CHROM_A", "START_A", "END_B", "NAME_A", "DUPGENES"]],
        on=["CHROM_A", "START_A", "END_B", "NAME_A"],
        how="left",
    )

    # in the case where there are multiple genes in one cell (HUGO1@ENSEMBL1; HUGO2@ENSEMBL2; ...)
    # split them in to separate comma-separated columns
    merged["GENE_A"] = merged["GENE_A"].map(split_multi)
    merged[["SYMBOL_A", "GENEID_A"]] = merged["GENE_A"].str.split(";", n=1, expand=True)
    merged["GENE_B"] = merged["GENE_B"].map(split_multi)
    merged[["SYMBOL_B", "GENEID_B"]] = merged["GENE_B"].str.split(";", n=1, expand=True)

    merged["DELGENES"] = merged["DELGENES"].map(split_multi)
    merged[["DEL_SYMBOLS", "DEL_GENEIDS"]] = merged["DELGENES"].str.split(
        ";", n=1, expand=True
    )
    merged["DUPGENES"] = merged["DUPGENES"].map(split_multi)
    merged[["DUP_SYMBOLS", "DUP_GENEIDS"]] = merged["DUPGENES"].str.split(
        ";", n=1, expand=True
    )

    # drop redundant columns
    merged = merged.drop(["GENE_A", "GENE_B", "DELGENES", "DUPGENES"], axis=1)

    return merged


def split_multi(s: str) -> str:
    """
    Split multi-gene annotations into separate symbol and ID components.

    Converts comma-separated gene annotations in format "SYMBOL1@ID1,SYMBOL2@ID2"
    into semicolon-separated format "SYMBOL1, SYMBOL2;ID1, ID2".

    :param s: String containing gene annotations or NaN/missing value
    :returns: Formatted string with symbols and IDs separated by semicolon
    """

    if pd.isna(s) or s == ".":
        return ".;."
    else:
        arr = s.split(",")
        symbols = ", ".join([elem.split("@")[0] for elem in arr])
        ids = ", ".join([elem.split("@")[1] for elem in arr])
        return symbols + ";" + ids


def filter_svs(
    df: pd.DataFrame,
    cosmic_fusion_gene_pairs_path: Path,
    onco_tsg_path: Path,
    sv_gnomad_cutoff: float = 0.001,
    large_sv_size: float = 1e9,
) -> pd.DataFrame:
    """
    Filter structural variants while rescuing clinically important ones.

    Applies size and frequency filters but rescues large SVs, variants affecting
    oncogenes/tumor suppressors, and variants creating known COSMIC fusion pairs.

    :param df: SVs in BEDPE format
    :param cosmic_fusion_gene_pairs_path: Path to file containing COSMIC fusion
        gene pairs
    :param onco_tsg_path: Path to file containing oncogene and tumor suppressor
        gene symbols
    :param sv_gnomad_cutoff: Maximum gnomAD allele frequency for an SV to be
        considered somatic
    :param large_sv_size: Size threshold beyond which SVs are considered large
        and need to be rescued
    :returns: Filtered SVs in BEDPE format
    """

    df["SVLEN_A"] = df["SVLEN_A"].astype("Int64")

    # drop variants shorter than 50
    df = df.loc[df["SVLEN_A"].isna() | df["SVLEN_A"].abs().ge(50)].copy()

    onco_tsg_df = pd.read_csv(
        onco_tsg_path, sep="\t", dtype="string", usecols=["hugo_symbol"]
    )
    oncogenes_and_ts = set(onco_tsg_df["hugo_symbol"])

    cosmic = pd.read_csv(cosmic_fusion_gene_pairs_path, dtype="string")
    cosmic_pairs = list(zip(cosmic["Gene_A"], cosmic["Gene_B"]))
    cosmic_pairs_sorted = set([tuple(sorted(elem)) for elem in cosmic_pairs])

    df["Rescue"] = False

    # rescue large SVs
    df.loc[df["SVLEN_A"].abs().ge(large_sv_size), "Rescue"] = True

    # rescue breakpoints that fall on oncogenes or tumor suppressors
    df["onco_ts_overlap_A"] = df["SYMBOL_A"].apply(
        onco_ts_overlap, oncogenes_and_ts=oncogenes_and_ts
    )
    df["onco_ts_overlap_B"] = df["SYMBOL_B"].apply(
        onco_ts_overlap, oncogenes_and_ts=oncogenes_and_ts
    )
    df.loc[df["onco_ts_overlap_A"] | df["onco_ts_overlap_A"], "Rescue"] = True

    # rescue gene pairs in cosmic
    df["pair_in_cosmic"] = df.apply(
        lambda row: list_all_pairs(
            row["SYMBOL_A"], row["SYMBOL_B"], cosmic_pairs_sorted
        ),
        axis=1,
    )
    df.loc[df["pair_in_cosmic"], "Rescue"] = True

    # gnomad AF parsing
    df["max_af"] = (
        df["vep_SV_overlap_AF_A"]
        .fillna("")
        .str.split("&")
        .apply(lambda x: max([float(e) if e != "" else 0 for e in x]))
    )

    # filter while keeping rescues
    df = df.loc[
        df["Rescue"] | df["max_af"].lt(sv_gnomad_cutoff),
        [
            "CHROM_A",
            "START_A",
            "END_A",
            "ID",
            "STRAND_A",
            "TYPE",
            "FILTER",
            "REF_A",
            "ALT_A",
            "SVLEN_A",
            "MATEID_A",
            "SVINSLEN_A",
            "BND_DEPTH_A",
            "MATE_BND_DEPTH_A",
            "SYMBOL_A",
            "GENEID_A",
            "vep_SV_overlap_name_A",
            "vep_SV_overlap_AF_A",
            "CHROM_B",
            "START_B",
            "END_B",
            "STRAND_B",
            "REF_B",
            "ALT_B",
            "SVLEN_B",
            "MATEID_B",
            "SVINSLEN_B",
            "BND_DEPTH_B",
            "MATE_BND_DEPTH_B",
            "SYMBOL_B",
            "GENEID_B",
            "vep_SV_overlap_name_B",
            "vep_SV_overlap_AF_B",
            "DEL_SYMBOLS",
            "DUP_SYMBOLS",
            "PR",
            "SR",
            "Rescue",
        ],
    ]

    return df


def onco_ts_overlap(s: str, oncogenes_and_ts: Set[str]) -> bool:
    """
    Check if any genes in a comma-separated list overlap with oncogenes/TSGs.

    :param s: Comma-separated string of gene symbols
    :param oncogenes_and_ts: Set of oncogene and tumor suppressor gene symbols
    :returns: True if any genes overlap with the oncogene/TSG set
    """

    l = s.split(", ")
    return len(set(l) & oncogenes_and_ts) > 0


def list_all_pairs(a: str, b: str, cosmic_pairs_sorted: Set[Tuple[str, str]]) -> bool:
    """
    Check if any gene pairs from two lists match known COSMIC fusion pairs.

    Creates all possible pairs from genes in lists a and b, then checks if any
    match the provided set of COSMIC fusion gene pairs.

    :param a: Comma-separated string of gene symbols from breakpoint A
    :param b: Comma-separated string of gene symbols from breakpoint B
    :param cosmic_pairs_sorted: Set of sorted tuples representing known COSMIC fusion
        gene pairs
    :returns: True if any gene pair matches a COSMIC fusion pair
    """

    alist = a.split(", ")
    blist = b.split(", ")

    all_pairs = list(itertools.product(alist, blist))
    all_pairs = set([tuple(sorted(elem)) for elem in all_pairs])

    return len(all_pairs & cosmic_pairs_sorted) > 0
