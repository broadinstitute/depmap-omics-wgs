#!/bin/zsh

set -euo pipefail

curl -SL https://snpeff.odsp.astrazeneca.com/databases/v5_2/snpEff_v5_2_GRCh38.mane.1.2.ensembl.zip \
  -o /tmp/snpEff_v5_2_GRCh38.mane.1.2.ensembl.zip

gcloud storage cp /tmp/snpEff_v5_2_GRCh38.mane.1.2.ensembl.zip gs://ccleparams/snpeff/
