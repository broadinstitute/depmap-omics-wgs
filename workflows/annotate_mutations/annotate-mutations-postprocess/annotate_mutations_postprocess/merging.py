import concurrent
import os
import re
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from csv import QUOTE_NONE
from glob import glob
from pathlib import Path

import bgzip
import pandas as pd
import tqdm
from click import echo


def do_merge(
    vcf_paths: list[Path], out_path: Path, chunk_size: int = 1000000000
) -> None:
    header_lines = get_header_lines(vcf_paths)

    with tempfile.TemporaryDirectory() as tmp_dir:
        chunks = make_chunks(header_lines, chunk_size)
        chunk_vcfs(vcf_paths, tmp_dir, chunks)
        process_chunks(header_lines, out_path, tmp_dir, chunks, len(vcf_paths))


def get_header_lines(vcf_paths: list[Path]) -> list[str]:
    header_lines = []
    col_header_line = None

    for path in vcf_paths:
        # don't read more lines than necessary to get the entire header
        break_next_time = False

        with open(path, "rb") as raw:
            echo(f"Reading {os.path.basename(path)} header")
            this_header_texts = ""  # start collecting header text

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
                # extract the line with column names
                col_header_line = [
                    x for x in this_header_lines if x.startswith("#CHROM")
                ][0]

            this_header_lines = [x for x in this_header_lines if x.startswith("##")]

            # add to the collected header lines
            header_lines.extend(this_header_lines)

    # de-dup but keep original order of lines
    return pd.Series([*header_lines, col_header_line]).drop_duplicates().tolist()


def make_chunks(header_lines: list[str], chunk_size: int) -> pd.DataFrame:
    # contig header lines contain names and lengths of chroms
    contig_lines = [x for x in header_lines if x.startswith("##contig")]
    chr_lengths_search = [
        re.search(r"^##contig=<ID=chr([^,]+),length=(\d+)>$", x) for x in contig_lines
    ]
    chr_lengths = pd.DataFrame(
        [x.groups() for x in chr_lengths_search if x is not None],
        columns=["chr", "length"],
    ).astype({"chr": "string", "length": "int64"})

    # remove things like decoy contigs
    valid_chrs = {str(x) for x in range(1, 23)}.union({"X", "Y"})
    chr_lengths = chr_lengths.loc[chr_lengths["chr"].isin(valid_chrs)]

    # start splitting regions into smaller chunks
    chunks = []

    for _, r in chr_lengths.iterrows():
        # `bcftools view` expects 1-index regions in TSVs used for splitting
        start = 1

        while start - 1 < r["length"]:
            end = min(start + chunk_size - 1, r["length"])
            chunks.append(
                {
                    "chr": f"chr{r["chr"]}",
                    "start": start,
                    "end": end,
                    "split_dir_name": f"chr{r['chr']}_{start}-{end}",
                }
            )
            start = end + 1

    return pd.DataFrame(chunks)


def chunk_vcfs(vcf_paths: list[Path], tmp_dir: str, chunks: pd.DataFrame) -> None:
    # need to make sure each VCF is tabix-indexed first
    with ProcessPoolExecutor() as executor:
        echo(f"Indexing {len(vcf_paths)} VCFs")
        futures = [executor.submit(index_vcf, path) for path in vcf_paths]

        for _ in tqdm.tqdm(
            concurrent.futures.as_completed(futures),
            total=len(vcf_paths),
            mininterval=1.0,
        ):
            pass

    # make the output split dirs
    chunks["split_dir"] = chunks["split_dir_name"].apply(
        lambda x: os.path.join(tmp_dir, x)
    )

    for _, r in chunks.iterrows():
        os.makedirs(r["split_dir"])

    # cross VCF paths with chunk paths
    vcf_chunks = pd.DataFrame({"path": vcf_paths}).merge(chunks, how="cross")
    vcf_chunks["chunk_path"] = vcf_chunks.apply(
        lambda x: os.path.join(x["split_dir"], os.path.basename(x["path"])), axis=1
    )
    vcf_chunks = vcf_chunks.drop(columns=["split_dir_name", "split_dir"])

    # write the chunks
    with ProcessPoolExecutor() as executor:
        echo(f"Writing {len(vcf_chunks)} VCF chunks")
        futures = [
            executor.submit(write_vcf_chunk, r)
            for r in vcf_chunks.to_dict(orient="records")
        ]

        for _ in tqdm.tqdm(
            concurrent.futures.as_completed(futures),
            total=len(vcf_chunks),
            mininterval=1.0,
        ):
            pass


def index_vcf(path: Path | str) -> None:
    subprocess.run(["bcftools", "index", path, "--force"])


def write_vcf_chunk(r: dict[str, str | int]) -> None:
    subprocess.run(
        [
            "bcftools",
            "view",
            r["path"],
            "--regions",
            f"{r['chr']}:{r['start']}-{r['end']}",
            "--output",
            r["chunk_path"],
        ]
    )


def process_chunks(
    header_lines: list[str],
    out_path: Path,
    tmp_dir: str,
    chunks: pd.DataFrame,
    n_vcfs: int,
) -> None:
    with open(out_path, "wb") as raw:
        with bgzip.BGZipWriter(raw) as f:
            # start writing bgzipped output VCF
            echo(f"Writing header to {out_path}")
            s = "\n".join(header_lines) + "\n"
            f.write(s.encode())

            for _, r in chunks.iterrows():
                split_dir = os.path.join(tmp_dir, r["split_dir_name"])

                # there should be one file per split per input VCF
                split_files = glob(
                    os.path.join(split_dir, "*.vcf.gz"), root_dir=tmp_dir
                )

                assert len(split_files) == n_vcfs, (
                    f"Expecting {n_vcfs} in {r['split_dir_name']}, "
                    f"but found {len(split_files)}"
                )

                echo(f"Merging {len(split_files)} VCF chunks for {r['split_dir_name']}")
                df = merge_chunks(split_files)

                echo(f"Writing merged {r['split_dir_name']} rows to {out_path}")
                s = df.to_csv(header=False, index=False, sep="\t", quoting=QUOTE_NONE)
                f.write(s.encode())


def merge_chunks(split_files: list[str]) -> pd.DataFrame:
    df = None

    # reduce VCF chunks
    with ThreadPoolExecutor() as executor:
        futures = [executor.submit(read_vcf, p) for p in split_files]

        for future in concurrent.futures.as_completed(futures):
            df = merge_two_vcfs(df, future.result())

    # de-dup each info field's annotations
    df["info"] = df["info"].apply(collect_uniq_annots).replace({"": "."})

    # df contains only one chrom, so sorting by pos is sufficient
    df = df.sort_values("pos")

    return df.fillna(".")


def read_vcf(path: str) -> pd.DataFrame:
    with open(path, "rb") as raw:
        with bgzip.BGZipReader(raw) as f:
            # noinspection PyTypeChecker
            df = pd.read_csv(
                f,
                sep="\t",
                header=None,
                names=[
                    "chrom",
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
                dtype="string",
                na_values=["."],
                keep_default_na=False,
                quoting=QUOTE_NONE,
                encoding_errors="backslashreplace",
            )

            # remove header rows
            df = df.loc[~df["chrom"].str.startswith("#")]
            df["pos"] = df["pos"].astype("int64")

            return df


def merge_two_vcfs(df1: pd.DataFrame | None, df2: pd.DataFrame) -> pd.DataFrame:
    if df1 is None:
        # starting base data frame
        df1 = df2
    else:
        # joining `id` and `info` fields to existing df
        df1 = df1.merge(
            df2,
            how="outer",
            on=["chrom", "pos", "ref", "alt"],
            suffixes=("", "2"),
        )

        # after correcting for variants potentially only being on one side of the merge,
        # ensure that columns that should be invariant are the same on each side
        for c in ["id", "qual", "filter", "format", "values"]:
            df1[c] = df1[c].fillna(df1[f"{c}2"])
            assert df1[c].eq(df1[f"{c}2"]).all()

        # concat info fields for each row (will de-dup later)
        df1["info"] = df1["info"].fillna("") + ";" + df1["info2"].fillna("")

        df1 = df1.drop(
            columns=["id2", "qual2", "filter2", "format2", "values2", "info2"]
        )

    return df1


def collect_uniq_annots(x: str) -> str:
    annots = {x for x in re.split(r";+", x) if len(x) > 0}
    annots = sorted(list(annots))
    return ";".join(annots)
