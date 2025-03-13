version 1.0

workflow postprocess_mutations {
    input {
        String workflow_version = "1.0"
        String workflow_source_url # populated automatically with URL of this script

        Array[File] duckdb
    }

    call postprocess {
        input:
            duckdb = duckdb
    }

    output {
        File somatic_variants = postprocess.somatic_variants
    }
}

task postprocess {
    input {
        Array[File] duckdb

        String docker_image
        String docker_image_hash_or_tag
        Int mem_gb = 12
        Int cpu = 2
        Int preemptible = 3
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
            --out-file="somatic_variants.parquet"
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
