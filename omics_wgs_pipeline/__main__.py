from pathlib import Path
from typing import Annotated, Optional

import pandas as pd
import tomllib
import typer

from omics_wgs_pipeline import terra
from omics_wgs_pipeline.terra_workflow import TerraWorkflow

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

config = {}


# noinspection PyUnusedLocal
def done(*args, **kwargs):
    typer.echo("Done.")


@app.callback(result_callback=done)
def main(config_path: Annotated[Optional[Path], typer.Option(exists=True)]):
    with open(config_path, "rb") as f:
        config.update(tomllib.load(f))


@app.command()
def update_workflow(workflow_name: Annotated[str, typer.Option()]) -> None:
    terra_workflow = TerraWorkflow(
        gcp_project_id=config["gcp_project_id"],
        pipelines_bucket_name=config["terra"]["pipelines_bucket_name"],
        repo_namespace=config["terra"]["repo_namespace"],
        repo_method_name=config["terra"][workflow_name]["repo_method_name"],
        workspace_namespace=config["terra"]["workspace_namespace"],
        workspace_name=config["terra"]["workspace_name"],
        method_config_name=config["terra"][workflow_name]["method_config_name"],
        method_synopsis=config["terra"][workflow_name]["method_synopsis"],
        workflow_wdl_path=Path(
            config["terra"][workflow_name]["workflow_wdl_path"]
        ).resolve(),
        method_config_json_path=Path(
            config["terra"][workflow_name]["method_config_json_path"]
        ).resolve(),
        firecloud_owners=config["terra"]["firecloud_owners"],
    )
    terra_workflow.update_workflow()


@app.command()
def refresh_terra_samples() -> None:
    terra.refresh_terra_samples(
        workspace_namespace=config["terra"]["workspace_namespace"],
        workspace_name=config["terra"]["workspace_name"],
    )


if __name__ == "__main__":
    app()
