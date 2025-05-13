version 1.0

workflow postprocess_mutations {
    input {
        String workflow_version = "1.0"
        String workflow_source_url # populated automatically with URL of this script

        Array[File] duckdb
        Float min_af = 0.15
        Int min_depth = 5
        Float max_pop_af = 0.00001
        Float max_brca1_func_assay_score = -1.328
    }

    call postprocess {
        input:
            duckdb = duckdb,
            min_af = min_af,
            min_depth = min_depth,
            max_pop_af = max_pop_af,
            max_brca1_func_assay_score = max_brca1_func_assay_score
    }

    output {
        File mut_somatic_variants = postprocess.somatic_variants
    }
}

task postprocess {
    input {
        Array[File] duckdb
        Float min_af
        Int min_depth
        Float max_pop_af
        Float max_brca1_func_assay_score

        String docker_image
        String docker_image_hash_or_tag
        Int mem_gb = 16
        Int cpu = 2
        Int preemptible = 2
        Int max_retries = 2
        Int additional_disk_gb = 0
    }

    Int disk_space = ceil(3 * size(duckdb, "GiB")) + 10 + additional_disk_gb

    command <<<
        set -euo pipefail

        mkdir -p parq
        mv ~{sep=" " duckdb} parq/

        python -m annotate_mutations_postprocess duckdb-to-maf \
            --db="tmp.duckdb" \
            --parquet-dir="./parq" \
            --out-file="somatic_variants.parquet" \
            --min-af=~{min_af} \
            --min-depth=~{min_depth} \
            --max-pop-af=~{max_pop_af} \
            --max-brca1-func-assay-score=~{max_brca1_func_assay_score}

    >>>

    output {
        File somatic_variants = "somatic_variants.parquet"
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
