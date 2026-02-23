import csv
import logging
from pathlib import Path
from typing import Annotated, Any

import pandas as pd
import typer

from compute_molecular_signatures.tasks import make_maf as do_make_maf

pd.set_option("display.max_columns", 30)
pd.set_option("display.max_colwidth", 50)
pd.set_option("display.max_info_columns", 30)
pd.set_option("display.max_info_rows", 20)
pd.set_option("display.max_rows", 20)
pd.set_option("display.max_seq_items", None)
pd.set_option("display.width", 200)
pd.set_option("expand_frame_repr", True)
pd.set_option("mode.chained_assignment", "warn")

app = typer.Typer(rich_markup_mode=None, pretty_exceptions_enable=False)

config: dict[str, Any] = {}


# noinspection PyUnusedLocal
def done(*args, **kwargs):
    logging.info("Done.")


@app.command()
def make_maf(
    muts: Annotated[
        Path,
        typer.Option(help="path the mut_sig Parquet file", exists=True),
    ],
    maf_out: Annotated[
        Path,
        typer.Option(help="path the output the mut_sig MAF"),
    ],
):
    logger = logging.getLogger()
    logger.addHandler(logging.StreamHandler())
    logger.setLevel(logging.INFO)

    maf = do_make_maf(muts_path=muts)
    maf.to_csv(maf_out, sep="\t", index=False, quoting=csv.QUOTE_NONE)


@app.callback(result_callback=done)
def main():
    logger = logging.getLogger()
    logger.addHandler(logging.StreamHandler())
    logger.setLevel(logging.INFO)


if __name__ == "__main__":
    app()
