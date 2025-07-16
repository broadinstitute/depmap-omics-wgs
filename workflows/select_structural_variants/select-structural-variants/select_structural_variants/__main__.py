import logging
import sys
from pathlib import Path
from typing import Annotated

import pandas as pd
import typer

from select_structural_variants.utils import do_select_structural_variants

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


def set_up_gcp_friendly_logging(level: int = logging.INFO) -> None:
    """
    Configure logging so that logs are routed to stdout/stderr based on severity,
    for compatibility with Google Cloud Logging, and are prepended by timestamps.

    :param level: log level to set
    """

    logger = logging.getLogger()
    logger.setLevel(level)
    logger.handlers.clear()

    # formatter for all log levels
    formatter = logging.Formatter(
        "%(asctime)s %(levelname)s - %(message)s", "%Y-%m-%d %H:%M:%S"
    )

    # handler for DEBUG and INFO → stdout
    stdout_handler = logging.StreamHandler(sys.stdout)
    stdout_handler.setLevel(logging.DEBUG)
    stdout_handler.addFilter(lambda x: x.levelno < logging.WARNING)
    stdout_handler.setFormatter(formatter)

    # handler for WARNING and above → stderr
    stderr_handler = logging.StreamHandler(sys.stderr)
    stderr_handler.setLevel(logging.WARNING)
    stderr_handler.setFormatter(formatter)

    logger.addHandler(stdout_handler)
    logger.addHandler(stderr_handler)


def main(
    input_bedpe: Annotated[Path, typer.Option(exists=True)],
    gene_annotation: Annotated[Path, typer.Option(exists=True)],
    del_annotation: Annotated[Path, typer.Option(exists=True)],
    dup_annotation: Annotated[Path, typer.Option(exists=True)],
    cosmic_fusion_gene_pairs: Annotated[Path, typer.Option(exists=True)],
    onco_tsg: Annotated[Path, typer.Option(exists=True)],
    out: Annotated[Path, typer.Option()],
) -> None:
    set_up_gcp_friendly_logging()

    do_select_structural_variants(
        input_bedpe_path=input_bedpe,
        gene_annotation_path=gene_annotation,
        del_annotation_path=del_annotation,
        dup_annotation_path=dup_annotation,
        cosmic_fusion_gene_pairs_path=cosmic_fusion_gene_pairs,
        onco_tsg_path=onco_tsg,
        out_path=out,
    )

    logging.info("Done.")


if __name__ == "__main__":
    typer.run(main)
