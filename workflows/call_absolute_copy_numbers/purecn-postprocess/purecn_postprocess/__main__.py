from pathlib import Path
from typing import Annotated

import pandas as pd
import typer
from click import echo

from purecn_postprocess.utils import collect_outputs

pd.set_option("display.max_columns", 30)
pd.set_option("display.max_colwidth", 50)
pd.set_option("display.max_info_columns", 30)
pd.set_option("display.max_info_rows", 20)
pd.set_option("display.max_rows", 20)
pd.set_option("display.max_seq_items", None)
pd.set_option("display.width", 200)
pd.set_option("expand_frame_repr", True)
pd.set_option("mode.chained_assignment", "warn")


def main(
    solution: Annotated[Path, typer.Option(exists=True)],
    loh: Annotated[Path, typer.Option(exists=True)],
    out: Annotated[Path, typer.Option()],
) -> None:
    collect_outputs(solution, loh, out)
    echo("Done.")


if __name__ == "__main__":
    typer.run(main)
