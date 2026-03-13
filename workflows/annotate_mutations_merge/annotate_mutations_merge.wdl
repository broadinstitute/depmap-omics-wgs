version 1.0

workflow annotate_mutations_merge {
    input {
        String sample_id
        Array[File] input_vcfs
    }

    call merge_info {
        input:
            vcfs = select_all(input_vcfs),
            output_file_base_name = sample_id + "_info_merged",
    }

    call index_vcf {
        input:
            vcf = merge_info.vcf_info_merged
    }

    call vcf_to_duckdb {
        input:
            vcf = merge_info.vcf_info_merged,
            vcf_index = index_vcf.vcf_index
    }

    output {
        File mut_annot_vcf = merge_info.vcf_info_merged
        File mut_annot_vcf_index = index_vcf.vcf_index
        Array[File] mut_duckdb = vcf_to_duckdb.duckdb
    }
}

task index_vcf {
    input {
        File vcf
        String index_format = "CSI"

        String docker_image = "us-central1-docker.pkg.dev/depmap-omics/terra-images/bcftools"
        String docker_image_hash_or_tag = ":production"
        Int mem_gb = 4
        Int cpu = 1
        Int preemptible = 2
        Int max_retries = 1
        Int additional_disk_gb = 0
    }

    Int disk_space = ceil(size(vcf, "GiB")) + 10 + additional_disk_gb

    String vcf_basename = basename(vcf)
    String index_format_option = if index_format == "CSI" then "--csi" else "--tbi"
    String index_file_ext = if index_format == "CSI" then ".csi" else ".tbi"

    command <<<
        set -euo pipefail

        bcftools index \
            "~{vcf}" \
            --output="~{vcf_basename}~{index_file_ext}" \
            ~{index_format_option}
    >>>

    output {
        File vcf_index = "~{vcf_basename}~{index_file_ext}"
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

task merge_info {
    input {
        Array[File] vcfs
        String output_file_base_name

        String docker_image = "us-central1-docker.pkg.dev/depmap-omics/terra-images/vcf-info-merger"
        String docker_image_hash_or_tag = ":production"
        Int mem_gb = 16
        Int cpu = 4
        Int preemptible = 1
        Int max_retries = 1
        Int additional_disk_gb = 0
    }

    Int disk_space = ceil(10 * size(vcfs, "GiB")) + 10 + additional_disk_gb

    command <<<
        python -m vcf_info_merger merge \
            --vcf ~{sep=" --vcf " vcfs} \
            --out="~{output_file_base_name}.vcf.gz"
    >>>

    output {
        File vcf_info_merged = "~{output_file_base_name}.vcf.gz"
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

task vcf_to_duckdb {
    input {
        File vcf
        File vcf_index
        Int batch_size = 1000000

        String docker_image = "us-central1-docker.pkg.dev/depmap-omics/terra-images/vcf-to-duckdb"
        String docker_image_hash_or_tag = ":production"
        Int mem_gb = 16
        Int cpu = 4
        Int preemptible = 0
        Int max_retries = 0
        Int additional_disk_gb = 0
    }

    Int disk_space = ceil(10 * size(vcf, "GiB")) + 30 + additional_disk_gb
    Int duck_db_memory_limit = floor(mem_gb * 0.75)

    command <<<
        set -euo pipefail

        mkdir -p "./data/tmp"

        python -m vcf_to_duckdb convert \
            --vcf="~{vcf}" \
            --tab="tmp.tsv" \
            --db="tmp.duckdb" \
            --parquet-dir="./parq" \
            --no-multiallelics \
            --config="/app/config.json" \
            --batch-size=~{batch_size} \
            --compression-level=15 \
            --tmp-dir="./data/tmp" \
            --memory-limit="~{duck_db_memory_limit}GB"
    >>>

    output {
        Array[File] duckdb = glob("parq/*.*")
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
