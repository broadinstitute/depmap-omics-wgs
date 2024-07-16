from _csv import QUOTE_NONE
from pathlib import Path

import bgzip
import pandas as pd
from click import echo
from pandas._testing import assert_frame_equal

from annotate_mutations_postprocess.vcf_utils import read_vcf


def do_merge(vcf_paths: list[Path], out_path: Path) -> None:
    header_lines = []
    col_header_line = None
    df = None
    common_col_names = None

    for path in vcf_paths:
        with open(path, "rb") as raw:
            echo(f"Reading {path}")
            this_header_texts = ""

            with bgzip.BGZipReader(raw) as f:
                while True:
                    if not (d := f.read(10 * 1024)):
                        break

                    text = d.tobytes().decode()
                    d.release()
                    this_header_texts += text

                    if "\n#CHROM" in this_header_texts:
                        break

            this_header_lines = this_header_texts.split("\n")

            if col_header_line is None:
                col_header_line = [
                    x for x in this_header_lines if x.startswith("#CHROM")
                ][0]

            this_header_lines = [x for x in this_header_lines if x.startswith("##")]
            header_lines.extend(this_header_lines)

            raw.seek(0)

            with bgzip.BGZipReader(raw) as f:
                this_df = read_vcf(f)

                if df is None:
                    df = this_df
                    common_col_names = df.columns.drop("info")
                else:
                    assert_frame_equal(df[common_col_names], this_df[common_col_names])
                    df["info"] = df["info"] + ";" + this_df["info"]

    header_lines = (
        pd.Series([*header_lines, col_header_line]).drop_duplicates().tolist()
    )

    echo(f"Writing to {out_path}")
    with open(out_path, "wb") as raw:
        with bgzip.BGZipWriter(raw) as f:
            s = "\n".join(
                [
                    *header_lines,
                    df.to_csv(header=False, index=False, sep="\t", quoting=QUOTE_NONE),
                ]
            )

            f.write(s.encode())
