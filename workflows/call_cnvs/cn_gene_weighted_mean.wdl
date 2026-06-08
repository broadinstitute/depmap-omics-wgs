version 1.0

workflow cn_gene_weighted_mean {
    input {
        String sample_id
        File read_cov_bin
        File protein_coding_genes_bed
    }

    call calc_cn_gene_weighted_mean {
        input:
            sample_id = sample_id,
            read_cov_bin = read_cov_bin,
            protein_coding_genes_bed = protein_coding_genes_bed
    }

    output {
        File cn_by_gene_weighted_mean = calc_cn_gene_weighted_mean.cn_by_gene_weighted_mean
    }
}

task calc_cn_gene_weighted_mean {
    input {
        String sample_id
        File read_cov_bin
        File protein_coding_genes_bed

        String docker_image = "us-central1-docker.pkg.dev/depmap-omics/terra-images/call_segments"
        String docker_image_hash_or_tag = ":production"
        Int mem_gb = 4
        Int cpu = 1
        Int preemptible = 2
        Int max_retries = 0
        Int additional_disk_gb = 0
    }

    Int disk_space = ceil(
        size(read_cov_bin, "GiB")
        + size(protein_coding_genes_bed, "GiB")
    ) + 10 + additional_disk_gb

    command <<<
        set -euo pipefail

        # calculate weighted-mean copy numbers for protein-coding gene
        # remove NA bins so they don't inflate the denominator
        sed '1d' "~{read_cov_bin}" \
            | awk '$11 != "NA"' \
            | awk 'BEGIN{OFS="\t"} { $11 = 2^$11; print }' \
            | bedtools intersect -a "~{protein_coding_genes_bed}" -b stdin -wao \
            | awk '{ print $0"\t"$17*$18 }' \
            | sort --key="1,1V" --key="2,2n" --key="3,3n" --key="4,4n" \
            | bedtools groupby -g 1,2,3,4 -c 18,19,19 -o sum,sum,count \
            | awk '{ if ($6 == 0) { print $0"\tNA" } else { print $0"\t"$6/($5+1) } }' \
            | awk 'BEGIN{OFS="\t"} { if ($8 == "NA") { print } else { $8 = log($8)/log(2); print } }' \
            > "~{sample_id}.cn_gene_weighted_mean.tsv"
    >>>

    output {
        File cn_by_gene_weighted_mean = "~{sample_id}.cn_gene_weighted_mean.tsv"
    }

    runtime {
        docker: "~{docker_image}~{docker_image_hash_or_tag}"
        memory: "~{mem_gb} GB"
        disks: "local-disk ~{disk_space} SSD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }

    meta {
        allowNestedInputs: true
    }
}
