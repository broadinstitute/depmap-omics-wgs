"""
Simple CLI to call SignatureAnalyzer, which is installed in editable mode in
`signature_analyzer.Dockerfile`
"""

import argparse

import pandas as pd
import signatureanalyzer as sa

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

# load the MAF
maf = pd.read_parquet(args.input_maf)

# get context counts for the sample
spectra_df = sa.spectra.get_spectra_from_maf(
    maf, hgfile=args.ref_2bit, reference=args.sa_ref, real_snps=True
)[1]

# load the reference signatures
ref_df, ref_idx = sa.utils.load_reference_signatures(args.sa_ref, verbose=False)
ref_df = ref_df.set_index("Somatic Mutation Type").iloc[:, :-2]

# run supervised NMF
res_supervised = sa.supervised_bnmf.supervised_ardnmf(
    spectra_df, ref_df, objective="poisson", max_iter=args.max_iter, verbose=True
)

# collect the wide data frame of results
sig_df = res_supervised["H"].T

# make long and drop summary values
sig_df = (
    sig_df.rename(columns={sig_df.columns[0]: "score"})
    .T.drop(columns=["max", "max_id", "max_norm"])
    .T.reset_index(names="sig_id")
    .rename_axis(None, axis=1)
)

sig_df.to_parquet(args.parquet_out)
