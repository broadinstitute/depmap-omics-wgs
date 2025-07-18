import json
import logging
import os
import tomllib
from pathlib import Path
from typing import Annotated, Any

import pandas as pd
import typer
from nebelung.terra_workspace import TerraWorkspace

import depmap_omics_wgs.data as data
from depmap_omics_wgs.types import GumboClient
from depmap_omics_wgs.utils import (
    get_hasura_creds,
    get_secret_from_sm,
    make_workflow_from_config,
)

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
    logger = logging.getLogger()
    logger.addHandler(logging.StreamHandler())
    logger.setLevel(logging.INFO)

    with open(config_path, "rb") as f:
        config.update(tomllib.load(f))

    def get_gumbo_client() -> GumboClient:
        if config["gumbo_env"] == "dev":
            return GumboClient(
                url="http://localhost:8080/v1/graphql",
                username="dogspa",
                headers={"X-Hasura-Admin-Secret": "secret"},
            )

        # get URL and password for Gumbo GraphQL API from secrets manager
        hasura_creds = get_hasura_creds(gumbo_env=config["gumbo_env"])

        return GumboClient(
            url=hasura_creds["url"],
            username="dogspa",
            headers={"X-Hasura-Admin-Secret": hasura_creds["password"]},
        )

    ctx.obj = {
        "terra_workspace": TerraWorkspace(
            workspace_namespace=config["terra"]["workspace_namespace"],
            workspace_name=config["terra"]["workspace_name"],
            owners=json.loads(os.environ["FIRECLOUD_OWNERS"]),
        ),
        "get_gumbo_client": get_gumbo_client,
    }


@app.command()
def update_workflow(
    ctx: typer.Context, workflow_name: Annotated[str, typer.Option()]
) -> None:
    # need a GitHub PAT for persisting WDL in gists
    github_pat = get_secret_from_sm(
        "projects/201811582504/secrets/github-pat-for-wdl-gists/versions/latest"
    )

    terra_workflow = make_workflow_from_config(
        config, workflow_name, github_pat=github_pat
    )

    ctx.obj["terra_workspace"].update_workflow(terra_workflow=terra_workflow)


@app.command()
def refresh_terra_samples(ctx: typer.Context) -> None:
    data.refresh_terra_samples(
        terra_workspace=ctx.obj["terra_workspace"],
        gumbo_client=ctx.obj["get_gumbo_client"](),
        ref_urls=config["ref"],
    )


@app.command()
def delta_job(
    ctx: typer.Context,
    workflow_name: Annotated[str, typer.Option()],
    entity_type: Annotated[str, typer.Option()],
    entity_set_type: Annotated[str, typer.Option()],
    entity_id_col: Annotated[str, typer.Option()],
    expression: Annotated[str, typer.Option()],
    input_col: Annotated[list[str] | None, typer.Option()] = None,
    output_col: Annotated[list[str] | None, typer.Option()] = None,
) -> None:
    input_cols_set = None
    output_cols_set = None

    if input_col is not None:
        input_cols_set = set(input_col)

    if output_col is not None:
        output_cols_set = set(output_col)

    ctx.obj["terra_workspace"].submit_delta_job(
        terra_workflow=make_workflow_from_config(config, workflow_name),
        entity_type=entity_type,
        entity_set_type=entity_set_type,
        entity_id_col=entity_id_col,
        expression=expression,
        dry_run=config["dry_run"],
        input_cols=input_cols_set,
        output_cols=output_cols_set,
        max_n_entities=1,
    )


@app.command()
def refresh_legacy_terra_samples(
    ctx: typer.Context,
    sample_set_id: Annotated[str, typer.Option()],
) -> None:
    data.refresh_legacy_terra_samples(
        terra_workspace=ctx.obj["terra_workspace"],
        legacy_terra_workspace=TerraWorkspace(
            workspace_namespace=config["terra"]["legacy_workspace_namespace"],
            workspace_name=config["terra"]["legacy_workspace_name"],
        ),
        sample_set_id=sample_set_id,
        gumbo_client=ctx.obj["get_gumbo_client"](),
    )


@app.command()
def onboard_aligned_bams(ctx: typer.Context) -> None:
    data.onboard_aligned_bams(
        terra_workspace=ctx.obj["terra_workspace"],
        gumbo_client=ctx.obj["get_gumbo_client"](),
        gcp_project_id=config["gcp_project_id"],
        dry_run=config["dry_run"],
    )


if __name__ == "__main__":
    app()
