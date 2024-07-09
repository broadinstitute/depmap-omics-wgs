from pathlib import Path
from typing import Annotated, List, Optional

import pandas as pd
import typer
from click import echo

from vcf_to_depmap.conversion import convert

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
    vcf: Annotated[Path, typer.Option(exists=True)],
    oncogenes_list: Annotated[Path, typer.Option(exists=True)],
    tsg_list: Annotated[Path, typer.Option(exists=True)],
    out: Annotated[Path, typer.Option()],
    drop_col: Annotated[Optional[List[str]], typer.Option()] = None,
    force_keep: Annotated[Optional[List[str]], typer.Option()] = None,
    compound_info_field: Annotated[Optional[List[str]], typer.Option()] = None,
    url_encoded_col_name_regex: Annotated[Optional[str], typer.Option()] = None,
) -> None:
    drop_cols = set(drop_col) if drop_col is not None else {}
    force_keeps = set(force_keep) if force_keep is not None else {}
    compound_info_fields = (
        set(compound_info_field) if compound_info_field is not None else {}
    )

    convert(
        vcf_path=vcf,
        oncogenes_list_path=oncogenes_list,
        tsg_list_path=tsg_list,
        out_path=out,
        drop_cols=drop_cols,
        force_keep=force_keeps,
        compound_info_fields=compound_info_fields,
        url_encoded_col_name_regex=url_encoded_col_name_regex,
    )

    echo("Done.")


if __name__ == "__main__":
    typer.run(main)
