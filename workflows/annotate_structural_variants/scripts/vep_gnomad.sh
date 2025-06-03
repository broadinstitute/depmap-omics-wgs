#!/bin/zsh

bcftools view \
  -f PASS \
  gs://gcp-public-data--gnomad/release/4.1/genome_sv/gnomad.v4.1.sv.sites.vcf.gz \
  -o ./data/gnomad.v4.1.sv.sites.pass.vcf.gz
bcftools index -t ./data/gnomad.v4.1.sv.sites.pass.vcf.gz
