"""
Download the most recent hg38 data set available from
https://cancer.sanger.ac.uk/cosmic/download/cancer-mutation-census
"""

import pandas as pd

VERSION = "99"

cmc = pd.read_csv(
    f"./data/cmc/CancerMutationCensus_AllData_v{VERSION}_GRCh38.tsv.gz",
    sep="\t",
    low_memory=False,
    usecols=[
        "Mutation genome position GRCh38",
        "GENOMIC_WT_ALLELE_SEQ",
        "GENOMIC_MUT_ALLELE_SEQ",
        "MUTATION_SIGNIFICANCE_TIER",
    ],
    dtype="string",
)

cmc = cmc.dropna(subset="Mutation genome position GRCh38")
cmc = cmc.loc[cmc["MUTATION_SIGNIFICANCE_TIER"].ne("Other")]
cmc["MUTATION_SIGNIFICANCE_TIER"] = cmc["MUTATION_SIGNIFICANCE_TIER"].astype("int8")

cmc[["chrom", "start_end"]] = cmc["Mutation genome position GRCh38"].str.split(
    ":", regex=False, n=1, expand=True
)
cmc["chrom"] = "chr" + cmc["chrom"]
cmc = cmc.drop(columns="Mutation genome position GRCh38")

cmc[["start", "end"]] = (
    cmc["start_end"].str.split("-", regex=False, n=1, expand=True).astype("int64")
)
cmc = cmc.drop(columns=["start_end", "end"])

cmc = cmc.rename(
    columns={
        "start": "pos",
        "GENOMIC_WT_ALLELE_SEQ": "ref",
        "GENOMIC_MUT_ALLELE_SEQ": "alt",
        "MUTATION_SIGNIFICANCE_TIER": "mutation_significance_tier",
    }
)

cmc = cmc[["chrom", "pos", "ref", "alt", "mutation_significance_tier"]]
cmc = cmc.dropna(subset=["chrom", "pos", "mutation_significance_tier"])

cmc_ok = cmc.loc[~cmc.isna().any(axis=1)].copy()
cmc_to_fill = cmc.loc[cmc.isna().any(axis=1)].copy()

missing = cmc_to_fill.copy()
missing["start"] = missing["pos"] - 2
missing["end"] = missing["pos"] - 1
missing = missing[["chrom", "start", "end"]].drop_duplicates()
missing.to_csv("./data/cmc/missing.bed", sep="\t", header=False, index=False)

"""
bedtools getfasta \
    -fi ./data/Homo_sapiens_assembly38.fasta \
    -bed ./data/cmc/missing.bed \
    -bedOut \
    > ./data/cmc/missing_left_alleles.bed
"""

missing_left_alleles = pd.read_csv(
    "./data/cmc/missing_left_alleles.bed",
    sep="\t",
    header=None,
    names=["chrom", "start", "end", "left_allele"],
)
missing_left_alleles["pos"] = missing_left_alleles["end"] + 1
missing_left_alleles = missing_left_alleles[["chrom", "pos", "left_allele"]]

cmc_to_fill = cmc_to_fill.merge(
    missing_left_alleles, how="outer", on=["chrom", "pos"], indicator=True
)

assert cmc_to_fill["_merge"].eq("both").all()

cmc_to_fill[["ref", "alt"]] = cmc_to_fill[["ref", "alt"]].fillna("")
cmc_to_fill["ref"] = cmc_to_fill["left_allele"] + cmc_to_fill["ref"]
cmc_to_fill["alt"] = cmc_to_fill["left_allele"] + cmc_to_fill["alt"]
cmc_to_fill["pos"] -= 1

cmc_to_fill = cmc_to_fill.drop(columns=["left_allele", "_merge"])

cmc_fixed = pd.concat([cmc_ok, cmc_to_fill])

cmc_fixed["chrom"] = cmc_fixed["chrom"].str.lstrip("chr")
cmc_fixed["chrom"] = cmc_fixed["chrom"].replace({"X": "23", "Y": "24"}).astype("int8")
cmc_fixed = cmc_fixed.sort_values(["chrom", "pos"])
cmc_fixed["chrom"] = cmc_fixed["chrom"].astype("string").replace({"23": "X", "24": "Y"})
cmc_fixed["chrom"] = "chr" + cmc_fixed["chrom"]

cmc_fixed.to_csv("./data/cmc/cosmic_cmc.tsv", sep="\t", index=False, header=False)

"""
echo '##INFO=<ID=CMC_TIER,Number=1,Type=String,Description="COSMIC CMC Mutation Significance Tier">' \
    > ./data/cmc/cosmic_cmc.hdr.vcf

bcftools convert \
    --columns=CHROM,POS,REF,ALT,- \
    --fasta-ref=./data/Homo_sapiens_assembly38.fasta \
    --tsv2vcf \
    ./data/cmc/cosmic_cmc.tsv \
    -o ./data/cmc/cosmic_cmc.vcf

bgzip ./data/cmc/cosmic_cmc.tsv -k -f
tabix ./data/cmc/cosmic_cmc.tsv.gz -s1 -b2 -e2 -f
bcftools annotate ./data/cmc/cosmic_cmc.vcf \
    --annotations=./data/cmc/cosmic_cmc.tsv.gz \
    --output=./data/cmc/cosmic_cmc_annot.vcf \
    --header-lines=./data/cmc/cosmic_cmc.hdr.vcf \
    --columns=CHROM,POS,REF,ALT,CMC_TIER

bcftools norm ./data/cmc/cosmic_cmc_annot.vcf \
    --check-ref=w \
    --fasta-ref=./data/Homo_sapiens_assembly38.fasta \
    --output=./data/cmc/cosmic_cmc_annot_norm.vcf.gz
bcftools index ./data/cmc/cosmic_cmc_annot_norm.vcf.gz

bcftools view -H ./data/cmc/cosmic_cmc_annot_norm.vcf.gz \
    | cut -f 1-2,4-5,8 \
    | sed -e 's/CMC_TIER=//' \
    > ./data/cmc/cosmic_cmc_v99.tsv

bgzip ./data/cmc/cosmic_cmc_v99.tsv -k -f
tabix ./data/cmc/cosmic_cmc_v99.tsv.gz -s1 -b2 -e2 -f

# e.g.
bcftools annotate input.vcf.gz \
    --annotations=./data/cmc/cosmic_cmc_v99.tsv.gz \
    --output=output.vcf.gz \
    --header-lines=./data/cmc/cosmic_cmc.hdr.vcf \
    --columns=CHROM,POS,REF,ALT,CMC_TIER
"""

# rename and upload final .vcf.gz and .tbi files
