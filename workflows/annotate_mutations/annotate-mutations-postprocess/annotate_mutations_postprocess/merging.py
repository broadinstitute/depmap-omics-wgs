import os
import tempfile
from csv import QUOTE_NONE
from glob import glob
from pathlib import Path

import bgzip
import numpy as np
import pandas as pd
from click import echo
from pandas._testing import assert_frame_equal


def do_merge(vcf_paths: list[Path], out_path: Path) -> None:
    header_lines = get_header_lines(vcf_paths)

    with tempfile.TemporaryDirectory() as tmp_dir:
        chunk_vcfs(vcf_paths, tmp_dir)
        process_chunks(header_lines, out_path, tmp_dir)


def get_header_lines(vcf_paths: list[Path]) -> list[str]:
    header_lines = []
    col_header_line = None

    for path in vcf_paths:
        break_next_time = False

        with open(path, "rb") as raw:
            echo(f"Reading {os.path.basename(path)} header")
            this_header_texts = ""

            # assume file is bgzipped
            with bgzip.BGZipReader(raw) as f:
                while True:
                    # read a small block of bytes at a time
                    if not (d := f.read(10 * 1024)):
                        break

                    # concat the latest chunk of text
                    text = d.tobytes().decode()
                    d.release()
                    this_header_texts += text

                    # check if we've reached the end of the header section and get one
                    # more chunk
                    if break_next_time:
                        break
                    elif "\n#CHROM" in this_header_texts:
                        break_next_time = True

            # extract the header lines and the column headers
            this_header_lines = this_header_texts.split("\n")

            if col_header_line is None:
                col_header_line = [
                    x for x in this_header_lines if x.startswith("#CHROM")
                ][0]

            this_header_lines = [x for x in this_header_lines if x.startswith("##")]

            # add to the collected header lines
            header_lines.extend(this_header_lines)

    return pd.Series([*header_lines, col_header_line]).drop_duplicates().tolist()


def chunk_vcfs(vcf_paths: list[Path], tmp_dir: str) -> None:
    for path in vcf_paths:
        with open(path, "rb") as raw:
            with bgzip.BGZipReader(raw) as f:
                echo(f"Reading {os.path.basename(path)}")
                # noinspection PyTypeChecker
                df = pd.read_csv(
                    f,
                    sep="\t",
                    header=None,
                    names=[
                        "chromosome",
                        "position",
                        "id",
                        "ref",
                        "alt",
                        "qual",
                        "filter",
                        "info",
                        "format",
                        "values",
                    ],
                    dtype="string",
                    keep_default_na=False,
                    quoting=QUOTE_NONE,
                    encoding_errors="backslashreplace",
                )

                # remove header rows
                df = df.loc[~df["chromosome"].str.startswith("#")]
                df["position"] = df["position"].astype("int64")

                # split into named chunks with sortable names
                df["chunk"] = 1 + np.arange(len(df)) // 100000
                n_chunks = df["chunk"].iloc[-1]
                df["chunk"] = df["chunk"].astype("string").str.zfill(len(str(n_chunks)))

                echo(f"Writing {os.path.basename(path)} as {n_chunks} parquet files")
                df.to_parquet(
                    path=tmp_dir,
                    partition_cols=["chunk"],
                    index=False,
                    engine="pyarrow",
                )


def process_chunks(header_lines: list[str], out_path: Path, tmp_dir: str) -> None:
    with open(out_path, "wb") as raw:
        with bgzip.BGZipWriter(raw) as f:
            echo(f"Writing header to {out_path}")
            s = "\n".join(header_lines) + "\n"
            f.write(s.encode())

            # get all split subfolders
            split_dirs = sorted(glob("*", root_dir=tmp_dir))

            for sub_dir in split_dirs:
                # there should be one Parquet file per subfolder
                parquet_files = glob(
                    os.path.join(sub_dir, "*.parquet"), root_dir=tmp_dir
                )

                df = None

                for p in parquet_files:
                    echo(f"Reading {p}")
                    this_df = pd.read_parquet(
                        os.path.join(tmp_dir, p), engine="pyarrow"
                    )

                    if df is None:
                        # starting base data frame
                        df = this_df
                    else:
                        # some annotators might populate the ID col, so collect those
                        df["id"] = df["id"].replace({".": pd.NA}).fillna(this_df["id"])

                        # most columns should match since annotators don't modify them
                        assert_frame_equal(
                            df[df.columns.drop(["id", "info"])],
                            this_df[this_df.columns.drop(["id", "info"])],
                        )

                        # concat info fields for each row
                        df["info"] = df["info"] + ";" + this_df["info"]

                # de-dup each info field's annotations
                df["info"] = df["info"].apply(
                    lambda x: ";".join(sorted(list(set(x.split(";")))))
                )

                s = df.to_csv(header=False, index=False, sep="\t", quoting=QUOTE_NONE)

                echo(f"Writing merged {sub_dir} rows to {out_path}")
                f.write(s.encode())
