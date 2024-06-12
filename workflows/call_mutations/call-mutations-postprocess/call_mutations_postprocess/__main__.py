from pathlib import Path
from typing import Annotated

import typer
from click import echo

from call_mutations_postprocess.utils import fix_clustered_event_flag


def main(
    in_vcf: Annotated[Path, typer.Option(exists=True)],
    out_vcf: Annotated[Path, typer.Option()],
    region_size: Annotated[int, typer.Option()],
) -> None:
    fix_clustered_event_flag(in_vcf, out_vcf, region_size)
    echo("Done.")


if __name__ == "__main__":
    typer.run(main)
