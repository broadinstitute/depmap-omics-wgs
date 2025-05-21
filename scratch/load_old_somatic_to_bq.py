from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd
from google.api_core.exceptions import NotFound
from google.cloud import bigquery
from nebelung.terra_workspace import TerraWorkspace

pd.set_option("display.max_columns", 30)
pd.set_option("display.max_colwidth", 50)
pd.set_option("display.max_info_columns", 30)
pd.set_option("display.max_info_rows", 20)
pd.set_option("display.max_rows", 20)
pd.set_option("display.max_seq_items", None)
pd.set_option("display.width", 200)
pd.set_option("expand_frame_repr", True)
pd.set_option("mode.chained_assignment", "warn")

bq_client = bigquery.Client(project="depmap-omics")

ws = TerraWorkspace(
    workspace_namespace="broad-firecloud-ccle",
    workspace_name="depmap-omics-wgs-mut-dev",
)

samples = ws.get_entities("sample")

eval_samples = samples.loc[
    :,
    ["sample_id", "wgs_cn_depmap_maf_25q2"],
].iloc[1:]

eval_samples.set_index("sample_id", inplace=True)


def upload_parquet_to_bq(table_id: str, gcs_uri: str):
    job_config = bigquery.LoadJobConfig(
        source_format=bigquery.SourceFormat.CSV,
        autodetect=True,
        write_disposition=bigquery.WriteDisposition.WRITE_APPEND,
        column_name_character_map="V2",
    )

    load_job = bq_client.load_table_from_uri(gcs_uri, table_id, job_config=job_config)
    load_job.result()
    print(f"Loaded {load_job.output_rows} rows into {table_id}")


with ThreadPoolExecutor(max_workers=4) as executor:
    futures = []

    for sample_id, f in eval_samples["wgs_cn_depmap_maf_25q2"].items():
        table_name = sample_id.replace("-", "_")
        table_id = f"{bq_client.project}.mut_eval.{table_name}_somatic"

        try:
            bq_client.get_table(table_id)
            print(f"Table {table_id} already exists. Skipping upload.")
            continue
        except NotFound:
            pass

        futures.append(
            executor.submit(upload_parquet_to_bq, table_id=table_id, gcs_uri=f)
        )

    for future in as_completed(futures):
        try:
            future.result()
        except Exception as e:
            print(f"Error: {e}")
