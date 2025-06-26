version 1.0

workflow select_somatic_variants {
    input {
        String sample_id
        Array[File] duckdb
        Float min_af = 0.15
        Int min_depth = 5
        Float max_pop_af = 0.00001
        Float max_brca1_func_assay_score = -1.328
    }

    call do_select_somatic_variants {
        input:
            sample_id = sample_id,
            duckdb = duckdb,
            min_af = min_af,
            min_depth = min_depth,
            max_pop_af = max_pop_af,
            max_brca1_func_assay_score = max_brca1_func_assay_score
    }

    output {
        File mut_enriched_variants = do_select_somatic_variants.enriched_variants
        File mut_somatic_variants = do_select_somatic_variants.somatic_variants
    }
}

task do_select_somatic_variants {
    input {
        String sample_id
        Array[File] duckdb
        Float min_af
        Int min_depth
        Float max_pop_af
        Float max_brca1_func_assay_score
        Int batch_size = 500000

        String docker_image
        String docker_image_hash_or_tag
        Int mem_gb = 32
        Int cpu = 8
        Int preemptible = 1
        Int max_retries = 0
        Int additional_disk_gb = 0
    }

    Int disk_space = ceil(10 * size(duckdb, "GiB")) + 20 + additional_disk_gb

    command <<<
        set -euo pipefail

        mkdir -p parq
        mv ~{sep=" " duckdb} parq/

        mkdir -p db_tmp

        python -m select_somatic_variants \
            --db="tmp.duckdb" \
            --parquet-dir="./parq" \
            --variants-enriched-out-file="~{sample_id}.enriched_variants.parquet" \
            --somatic-variants-out-file="~{sample_id}.somatic_variants.parquet" \
            --db-tmp-dir-path="./db_tmp/" \
            --sample-id="~{sample_id}" \
            --min-af=~{min_af} \
            --min-depth=~{min_depth} \
            --max-pop-af=~{max_pop_af} \
            --max-brca1-func-assay-score=~{max_brca1_func_assay_score} \
            --batch-size=~{batch_size}
    >>>

    output {
        File enriched_variants = "~{sample_id}.enriched_variants.parquet"
        File somatic_variants = "~{sample_id}.somatic_variants.parquet"
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
