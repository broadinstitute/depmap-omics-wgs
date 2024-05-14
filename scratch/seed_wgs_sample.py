import pandas as pd

pd.set_option("display.max_columns", 30)
pd.set_option("display.max_colwidth", 50)
pd.set_option("display.max_info_columns", 30)
pd.set_option("display.max_info_rows", 20)
pd.set_option("display.max_rows", 20)
pd.set_option("display.max_seq_items", None)
pd.set_option("display.width", 200)
pd.set_option("expand_frame_repr", True)
pd.set_option("mode.chained_assignment", "warn")

samples = pd.read_table("./data/DepMap_WGS_CN-sample.tsv", sep="\t", low_memory=False)
samples = samples.convert_dtypes()

cram_bam_cols = [
    col
    for col in samples.columns
    if samples[col].astype("string").str.contains(r"\.(cram|bam)$", regex=True).any()
]

cram_bams = samples.loc[:, ["entity:sample_id", *cram_bam_cols]]

samples = samples.rename(
    columns={
        "entity:sample_id": "sequencing_id",
        "hg38_cram_filepath": "delivered_bam_cram",
        "hg38_crai_filepath": "delivered_bai_crai",
        "hg19_bam_filepath": "delivered_bam_cram_hg19",
        "hg19_bai_filepath": "delivered_bai_crai_hg19",
        "duplication_metrics": "mark_dup_metrics",
        "internal_bam_filepath": "bam",
        "internal_bai_filepath": "bai",
    }
)

samples = samples.loc[
    :,
    [
        "sequencing_id",
        "delivered_bam_cram",
        "delivered_bai_crai",
        "delivered_bam_cram_hg19",
        "delivered_bai_crai_hg19",
        "mark_dup_metrics",
        "bam",
        "bai",
    ],
]

samples["delivered_bam_cram"] = samples["delivered_bam_cram"].fillna(
    samples["delivered_bam_cram_hg19"]
)
samples["delivered_bai_crai"] = samples["delivered_bai_crai"].fillna(
    samples["delivered_bai_crai_hg19"]
)

samples["delivered_reference_genome_id"] = "hg38"
samples.loc[
    samples["delivered_bam_cram"].eq(samples["delivered_bam_cram_hg19"]),
    "delivered_reference_genome_id",
] = "hg19"

samples = samples.drop(columns=["delivered_bam_cram_hg19", "delivered_bai_crai_hg19"])

samples = samples.dropna(subset=["bam", "bai"])

samples["created_at"] = pd.Timestamp.now(tz="UTC").isoformat()
samples["updated_at"] = samples["created_at"]
samples["reference_genome_id"] = "hg38"

samples = samples.loc[~samples["sequencing_id"].isin({"CDS-fUSLRB"})]

samples.to_csv("./data/wgs_samples.csv", index=False)
