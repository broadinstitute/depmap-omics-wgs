import logging
from pathlib import Path
from typing import Annotated, Any, List

import pandas as pd
import tomllib
import typer
from vcf_info_merger import info_merge_vcfs

from annotate_mutations_postprocess.annotation import annotate_vcf

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

config: dict[str, Any] = {}


# noinspection PyUnusedLocal
def done(*args, **kwargs):
    logging.info("Done.")


@app.callback(result_callback=done)
def main(
    ctx: typer.Context,
    config_path: Annotated[Path, typer.Option(exists=True)],
):
    logging.basicConfig(
        format="%(asctime)s %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=logging.INFO,
    )

    with open(config_path, "rb") as f:
        config.update(tomllib.load(f))

    ctx.obj = config


@app.command()
def merge_info(
    vcf: Annotated[List[Path], typer.Option(exists=True)],
    out: Annotated[Path, typer.Option()],
) -> None:
    info_merge_vcfs(vcf_paths=vcf, out_path=out)


@app.command()
def vcf_to_depmap(
    db: Annotated[Path, typer.Option()],
    parquet_dir: Annotated[Path, typer.Option()],
    oncogenes: Annotated[Path, typer.Option(exists=True)],
    tumor_suppressor_genes: Annotated[Path, typer.Option(exists=True)],
) -> None:
    with open(oncogenes, "r") as f:
        oncogenes_list = set(line.strip() for line in f)

    with open(tumor_suppressor_genes, "r") as f:
        tumor_suppressor_genes_list = set(line.strip() for line in f)

    annotate_vcf(
        db_path=db,
        parquet_dir_path=parquet_dir,
        oncogenes=oncogenes_list,
        tumor_suppressor_genes=tumor_suppressor_genes_list,
    )


if __name__ == "__main__":
    app()
