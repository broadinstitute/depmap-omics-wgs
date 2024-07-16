import re
from pathlib import Path
from typing import Annotated, Any, List

import pandas as pd
import tomllib
import typer
from click import echo

from annotate_mutations_postprocess.conversion import convert
from annotate_mutations_postprocess.merging import do_merge

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
    typer.echo("Done.")


@app.callback(result_callback=done)
def main(
    ctx: typer.Context,
    config_path: Annotated[Path, typer.Option(exists=True)],
):
    with open(config_path.name, "rb") as f:
        config.update(tomllib.load(f))

    ctx.obj = config


@app.command()
def merge_info(
    vcf: Annotated[List[Path], typer.Option(exists=True)],
    out: Annotated[Path, typer.Option()],
) -> None:
    do_merge(vcf_paths=vcf, out_path=out)


@app.command()
def clean_and_annotate(
    ctx: typer.Context,
    vcf: Annotated[Path, typer.Option(exists=True)],
    dna_repair_genes: Annotated[Path, typer.Option(exists=True)],
    oncogenes_list: Annotated[Path, typer.Option(exists=True)],
    tsg_list: Annotated[Path, typer.Option(exists=True)],
    out: Annotated[Path, typer.Option()],
) -> None:
    convert(
        vcf_path=vcf,
        dna_repair_genes_path=dna_repair_genes,
        oncogenes_list_path=oncogenes_list,
        tsg_list_path=tsg_list,
        out_path=out,
        drop_cols=set(ctx.obj["settings"]["drop_cols"]),
        na_cols=set(ctx.obj["settings"]["na_cols"]),
        bool_cols=set(ctx.obj["settings"]["bool_cols"]),
        force_keep=set(ctx.obj["settings"]["force_keep"]),
        compound_info_fields=set(ctx.obj["settings"]["compound_info_fields"]),
        url_encoded_col_name_regex=re.compile(
            "|".join(ctx.obj["settings"]["url_encoded_col_name_regexes"])
        ),
        funco_sanitized_col_name_regex=re.compile(
            "|".join(ctx.obj["settings"]["funco_sanitized_col_name_regexes"])
        ),
    )

    echo("Done.")


if __name__ == "__main__":
    app()
