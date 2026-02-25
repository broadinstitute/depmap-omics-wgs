"""
Simple CLI to call SignatureAnalyzer, which is installed in editable mode in
`signature_analyzer.Dockerfile`
"""

import argparse
import logging
import sys
import warnings

import pandas as pd

# Suppress SyntaxWarnings from twobitreader dependency during import
with warnings.catch_warnings():
    warnings.simplefilter("ignore", SyntaxWarning)
    import signatureanalyzer as sa


def set_up_gcp_friendly_logging(level: int = logging.INFO) -> None:
    """
    Configure logging so that logs are routed to stdout/stderr based on severity,
    for compatibility with Google Cloud Logging, and are prepended by timestamps.

    :param level: log level to set
    """

    logger = logging.getLogger()
    logger.setLevel(level)
    logger.handlers.clear()

    # formatter for all log levels
    formatter = logging.Formatter(
        "%(asctime)s %(levelname)s - %(message)s", "%Y-%m-%d %H:%M:%S"
    )

    # handler for DEBUG and INFO → stdout
    stdout_handler = logging.StreamHandler(sys.stdout)
    stdout_handler.setLevel(logging.DEBUG)
    stdout_handler.addFilter(lambda x: x.levelno < logging.WARNING)
    stdout_handler.setFormatter(formatter)

    # handler for WARNING and above → stderr
    stderr_handler = logging.StreamHandler(sys.stderr)
    stderr_handler.setLevel(logging.WARNING)
    stderr_handler.setFormatter(formatter)

    logger.addHandler(stdout_handler)
    logger.addHandler(stderr_handler)


set_up_gcp_friendly_logging()

parser = argparse.ArgumentParser()

parser.add_argument(
    "--input-maf", type=str, required=True, help="MAF file of the variants"
)

parser.add_argument(
    "--ref-2bit", type=str, required=True, help="2bit file for the reference genome"
)

parser.add_argument(
    "--sa-ref",
    type=str,
    required=True,
    help="SignatureAnalyzer reference signature to use",
)

parser.add_argument(
    "--max-iter",
    type=int,
    default=30000,
    required=False,
    help="Max iterations to run",
)

parser.add_argument(
    "--parquet-out", type=str, required=True, help="path to write output Parquet file"
)

args = parser.parse_args()

logging.info("Loading the MAF")
maf = pd.read_parquet(args.input_maf)

logging.info("Getting context counts for the sample")
spectra_df = sa.spectra.get_spectra_from_maf(
    maf, hgfile=args.ref_2bit, reference=args.sa_ref, real_snps=True
)[1]

logging.info("Loading the reference signatures")
ref_df, ref_idx = sa.utils.load_reference_signatures(args.sa_ref, verbose=False)
ref_df = ref_df.set_index("Somatic Mutation Type").iloc[:, :-2]

logging.info("Running supervised NMF")
res_supervised = sa.supervised_bnmf.supervised_ardnmf(
    spectra_df, ref_df, objective="poisson", max_iter=args.max_iter, verbose=True
)

logging.info("Collecting the wide data frame of results")
sig_df = res_supervised["H"].T

logging.info("Making long and dropping summary values")
sig_df = (
    sig_df.rename(columns={sig_df.columns[0]: "score"})
    .T.drop(columns=["max", "max_id", "max_norm"])
    .T.reset_index(names="sig_id")
    .rename_axis(None, axis=1)
)

logging.info("Writing output to parquet file")
sig_df.to_parquet(args.parquet_out)
