from pathlib import Path
from typing import Annotated

import typer
from click import echo

from purecn_postprocess.utils import collect_outputs


def main(
    solution: Annotated[Path, typer.Option(exists=True)],
    loh: Annotated[Path, typer.Option(exists=True)],
    out: Annotated[Path, typer.Option()],
) -> None:
    collect_outputs(solution, loh, out)
    echo("Done.")


if __name__ == "__main__":
    typer.run(main)
