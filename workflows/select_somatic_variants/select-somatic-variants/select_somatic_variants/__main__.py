import logging
import sys
from pathlib import Path
from typing import Annotated

import pandas as pd
import typer

from select_somatic_variants.utils import do_select_somatic_variants

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
    sample_id: Annotated[str, typer.Option()],
    db: Annotated[Path, typer.Option()],
    parquet_dir: Annotated[Path, typer.Option()],
    somatic_variants_out_file: Annotated[Path, typer.Option()],
    variants_enriched_out_file: Annotated[Path | None, typer.Option()] = None,
    db_tmp_dir_path: Annotated[Path, typer.Option()] = Path("/tmp/duckdb"),
    min_af: Annotated[float, typer.Option()] = 0.15,
    min_depth: Annotated[int, typer.Option()] = 5,
    max_pop_af: Annotated[float, typer.Option()] = 1e-05,
    max_brca1_func_assay_score: Annotated[float, typer.Option()] = -1.328,
    batch_size: Annotated[int, typer.Option()] = 1000000,
) -> None:
    set_up_gcp_friendly_logging()
    
    do_select_somatic_variants(
        db_path=db,
        parquet_dir_path=parquet_dir,
        variants_enriched_out_file_path=variants_enriched_out_file,
        somatic_variants_out_file_path=somatic_variants_out_file,
        db_tmp_dir_path=db_tmp_dir_path,
        sample_id=sample_id,
        min_af=min_af,
        min_depth=min_depth,
        max_pop_af=max_pop_af,
        max_brca1_func_assay_score=max_brca1_func_assay_score,
        batch_size=batch_size,
    )
    
    done()


if __name__ == "__main__":
    typer.run(main)
