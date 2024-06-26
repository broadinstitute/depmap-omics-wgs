import pandas as pd

hess = pd.read_csv(
    "./data/hess_drivers.csv",
    # index_col=True,
    usecols=[
        "CHROM",
        "POS",
        "REF",
        "ALT",
        "sig_assignments",
    ],
)

hess = hess.rename(columns={"sig_assignments": "HESS"})
hess["HESS"] = hess["HESS"].fillna("")
hess = hess.dropna()

hess = hess[["CHROM", "POS", "REF", "ALT", "HESS"]]

hess["chr_num"] = hess["CHROM"].str.lstrip("chr")
hess.loc[hess["chr_num"].eq("X"), "chr_num"] = "23"
hess.loc[hess["chr_num"].eq("Y"), "chr_num"] = "24"
hess["chr_num"] = hess["chr_num"].astype("int8")

hess = hess.sort_values(["chr_num", "POS"])
hess = hess.drop(columns="chr_num")

hess.to_csv("./data/hess_drivers.tsv", sep="\t", index=False, header=False)

"""
echo '##INFO=<ID=HESS,Number=1,Type=String,Description="Hess driver signature">' \
    > ./data/hess.hdr.vcf
bgzip ./data/hess_drivers.tsv -k
tabix ./data/hess_drivers.tsv.gz -s1 -b2 -e2
bcftools annotate input.vcf.gz \
    --annotations=./data/hess_drivers.tsv.gz \
    --output=./data/output.vcf.gz \
    --header-lines=./data/hess.hdr.vcf \
    --columns=CHROM,POS,REF,ALT,HESS
"""
