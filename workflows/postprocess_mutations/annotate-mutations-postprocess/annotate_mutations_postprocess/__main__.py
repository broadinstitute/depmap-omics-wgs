import logging
from pathlib import Path
from typing import Annotated

import pandas as pd
import typer

from annotate_mutations_postprocess.maf import convert_duckdb_to_maf

pd.set_option("display.max_columns", 30)
pd.set_option("display.max_colwidth", 50)
pd.set_option("display.max_info_columns", 30)
pd.set_option("display.max_info_rows", 20)
pd.set_option("display.max_rows", 20)
pd.set_option("display.max_seq_items", None)
pd.set_option("display.width", 200)
pd.set_option("expand_frame_repr", True)
pd.set_option("mode.chained_assignment", "warn")

app = typer.Typer()


# noinspection PyUnusedLocal
def done(*args, **kwargs):
    logging.info("Done.")


@app.callback(result_callback=done)
def main():
    logging.basicConfig(
        format="%(asctime)s %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=logging.INFO,
    )


@app.command()
def duckdb_to_maf(
    db: Annotated[Path, typer.Option()],
    parquet_dir: Annotated[Path, typer.Option()],
    out_file: Annotated[Path, typer.Option()],
    min_af: Annotated[float, typer.Option()] = 0.15,
    min_depth: Annotated[int, typer.Option()] = 2,
    max_pop_af: Annotated[float, typer.Option()] = 1e-05,
    max_brca1_func_assay_score: Annotated[float, typer.Option()] = -1.328,
) -> None:
    convert_duckdb_to_maf(
        db_path=db,
        parquet_dir_path=parquet_dir,
        out_file_path=out_file,
        min_af=min_af,
        min_depth=min_depth,
        max_pop_af=max_pop_af,
        max_brca1_func_assay_score=max_brca1_func_assay_score,
    )


if __name__ == "__main__":
    app()
