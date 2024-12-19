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
) -> None:
    convert_duckdb_to_maf(
        db_path=db, parquet_dir_path=parquet_dir, out_file_path=out_file
    )


if __name__ == "__main__":
    app()
