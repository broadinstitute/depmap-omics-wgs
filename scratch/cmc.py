import pandas as pd

cmc = pd.read_csv(
    "./data/CancerMutationCensus_AllData_v99_GRCh38.tsv.gz",
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

cmc = cmc.reset_index(drop=True)
cmc["id"] = cmc.index.values
cmc["id"].astype("string")

cmc = cmc[["chrom", "pos", "id", "ref", "alt", "mutation_significance_tier"]]

cmc["chrom"] = cmc["chrom"].replace({"X": "23", "Y": "24"}).astype("int8")
cmc = cmc.sort_values(["chrom", "pos"])
cmc["chrom"] = cmc["chrom"].astype("string").replace({"23": "X", "24": "Y"})
cmc["chrom"] = "chr" + cmc["chrom"]

cmc = cmc.dropna(subset=["chrom", "pos", "ref", "mutation_significance_tier"])
cmc["alt"] = cmc["alt"].fillna(".")

cmc.to_csv("./data/cmc/cosmic_cmc.tsv", sep="\t", index=False, header=False)

"""
echo '##INFO=<ID=CMC_TIER,Number=1,Type=String,Description="COSMIC CMC Mutation Significance Tier">' \
    > ./data/cmc/cosmic_cmc.hdr.vcf

bcftools convert \
    --columns=CHROM,POS,REF,ALT,- \
    --fasta-ref=./data//Homo_sapiens_assembly38.fasta \
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

gatk LeftAlignAndTrimVariants \
    -R ./data/Homo_sapiens_assembly38.fasta \
    -V ./data/cmc/cosmic_cmc_annot.vcf \
    -O ./data/cmc/cosmic_cmc_annot_norm.vcf.gz \
    --split-multi-allelics \
    --dont-trim-alleles
    
# e.g.
bcftools annotate input.vcf.gz \
    --annotations=./data/cmc/cosmic_cmc_annot_norm.vcf.gz \
    --output=output.vcf.gz \
    --header-lines=./data/cmc/cosmic_cmc.hdr.vcf \
    --columns=CMC_TIER
"""
