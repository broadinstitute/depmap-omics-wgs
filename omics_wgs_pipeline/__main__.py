import os
from pathlib import Path
from typing import Annotated, Any

import pandas as pd
import tomllib
import typer

from omics_wgs_pipeline.data import (
    delta_preprocess_wgs_samples,
    make_terra_samples,
    put_task_results,
)
from omics_wgs_pipeline.terra import TerraWorkflow, TerraWorkspace
from omics_wgs_pipeline.types import GumboClient

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

    ctx.obj = {
        "terra_workspace": TerraWorkspace(
            workspace_namespace=config["terra"]["workspace_namespace"],
            workspace_name=config["terra"]["workspace_name"],
            firecloud_owners=config["terra"]["firecloud_owners"],
        ),
        "gumbo_client": GumboClient(
            url=os.environ["HASURA_URL"],
            username="omics_wgs_pipeline",
            headers={"X-Hasura-Admin-Secret": os.environ["HASURA_ADMIN_SECRET"]},
        ),
    }


@app.command()
def update_workflow(
    ctx: typer.Context, workflow_name: Annotated[str, typer.Option()]
) -> None:
    terra_workflow = TerraWorkflow(
        gcp_project_id=config["gcp_project_id"],
        pipelines_bucket_name=config["terra"]["pipelines_bucket_name"],
        repo_namespace=config["terra"]["repo_namespace"],
        repo_method_name=config["terra"][workflow_name]["repo_method_name"],
        method_config_name=config["terra"][workflow_name]["method_config_name"],
        method_synopsis=config["terra"][workflow_name]["method_synopsis"],
        workflow_wdl_path=Path(
            config["terra"][workflow_name]["workflow_wdl_path"]
        ).resolve(),
        method_config_json_path=Path(
            config["terra"][workflow_name]["method_config_json_path"]
        ).resolve(),
    )

    ctx.obj["terra_workspace"].update_workflow(terra_workflow)


@app.command()
def refresh_terra_samples(ctx: typer.Context) -> None:
    samples = make_terra_samples(ctx.obj["gumbo_client"])
    ctx.obj["terra_workspace"].upload_entities(samples)


@app.command()
def run_workflow(
    ctx: typer.Context, workflow_name: Annotated[str, typer.Option()]
) -> None:
    terra_workflow = TerraWorkflow(
        gcp_project_id=config["gcp_project_id"],
        pipelines_bucket_name=config["terra"]["pipelines_bucket_name"],
        repo_namespace=config["terra"]["repo_namespace"],
        repo_method_name=config["terra"][workflow_name]["repo_method_name"],
        method_config_name=config["terra"][workflow_name]["method_config_name"],
        method_synopsis=config["terra"][workflow_name]["method_synopsis"],
        workflow_wdl_path=Path(
            config["terra"][workflow_name]["workflow_wdl_path"]
        ).resolve(),
        method_config_json_path=Path(
            config["terra"][workflow_name]["method_config_json_path"]
        ).resolve(),
    )

    if workflow_name == "preprocess_wgs_sample":
        delta_preprocess_wgs_samples(ctx.obj["terra_workspace"], terra_workflow)
    else:
        raise NotImplementedError(f"Workflow {workflow_name} not implemented")


@app.command()
def persist_outputs_in_gumbo(ctx: typer.Context) -> None:
    outputs = ctx.obj["terra_workspace"].collect_workflow_outputs()
    put_task_results(ctx.obj["gumbo_client"], config["gcp_project_id"], outputs)


if __name__ == "__main__":
    app()
