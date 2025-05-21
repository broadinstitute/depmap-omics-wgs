from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd
from google.cloud import bigquery, storage
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
    ["sample_id", "mut_enriched_variants"],
].dropna()

eval_samples.set_index("sample_id", inplace=True)


def upload_parquet_to_bq(table_id: str, gcs_uri: str):
    job_config = bigquery.LoadJobConfig(
        source_format=bigquery.SourceFormat.PARQUET,
        autodetect=True,
        write_disposition=bigquery.WriteDisposition.WRITE_APPEND,
        clustering_fields=["sample_id", "somatic", "chrom"],
    )

    load_job = bq_client.load_table_from_uri(gcs_uri, table_id, job_config=job_config)
    load_job.result()
    print(f"Loaded {load_job.output_rows} rows into {table_id}")


with ThreadPoolExecutor(max_workers=4) as executor:
    futures = [
        executor.submit(
            upload_parquet_to_bq,
            table_id=f"{bq_client.project}.mut_eval.new",
            gcs_uri=f,
        )
        for sample_id, f in eval_samples["mut_enriched_variants"].items()
    ]

    for future in as_completed(futures):
        try:
            future.result()
        except Exception as e:
            print(f"Error: {e}")
