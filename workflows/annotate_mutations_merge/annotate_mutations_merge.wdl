version 1.0

workflow annotate_mutations_merge {
    input {
        String sample_id
        Array[File] input_vcfs
        File xy_intervals
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

    call split_vcf_by_chrom {
        input:
            vcf = merge_info.vcf_info_merged,
            xy_intervals = xy_intervals
    }

    Array[Pair[File, File]] indexed_vcfs = zip(split_vcf_by_chrom.vcfs, split_vcf_by_chrom.indexes)

    scatter (indexed_vcf in indexed_vcfs) {
        call vcf_to_duckdb {
            input:
                vcf = indexed_vcf.left,
                vcf_index = indexed_vcf.right
        }
    }

    call merge_dbs {
        input:
            duckdbs = vcf_to_duckdb.duckdb
    }

    output {
        File mut_annot_vcf = merge_info.vcf_info_merged
        File mut_annot_vcf_index = index_vcf.vcf_index
        Array[File] mut_duckdb = merge_dbs.duckdb
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

task split_vcf_by_chrom {
    input {
        File vcf
        File xy_intervals

        String docker_image = "us-central1-docker.pkg.dev/depmap-omics/terra-images/bcftools"
        String docker_image_hash_or_tag = ":production"
        Int mem_gb = 2
        Int cpu = 1
        Int preemptible = 2
        Int max_retries = 1
        Int additional_disk_gb = 0
    }

    Int disk_space = ceil(2 * size(vcf, "GiB")) + 10 + additional_disk_gb

    command <<<
        set -euo pipefail

        bcftools index ~{vcf}

        mkdir "out"

        for chr in $(cat ~{xy_intervals}); do
            vcf_out="out/${chr}.vcf.gz"
            index_out="out/${chr}.vcf.gz.csi"

            bcftools view \
                "~{vcf}" \
                --regions="${chr}" \
                --output="${vcf_out}" \
                --no-version

            if [[ $(bcftools view -H "${vcf_out}" | wc -l) -eq 0 ]]; then
                rm "${vcf_out}"
            else
                bcftools index "${vcf_out}" -o "${index_out}"
            fi
        done
    >>>

    output {
        Array[File] vcfs = glob("out/*.vcf.gz")
        Array[File] indexes = glob("out/*.vcf.gz.csi")
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
        Int batch_size = 100000

        String docker_image = "us-central1-docker.pkg.dev/depmap-omics/terra-images/vcf-to-duckdb"
        String docker_image_hash_or_tag = ":production"
        Int mem_gb = 8
        Int cpu = 2
        Int preemptible = 1
        Int max_retries = 1
        Int additional_disk_gb = 0
    }

    Int disk_space = ceil(3 * size(vcf, "GiB")) + 10 + additional_disk_gb

    command <<<
        set -euo pipefail

        mkdir -p db_tmp

        python -m vcf_to_duckdb convert \
            --vcf="~{vcf}" \
            --tab="tmp.tsv" \
            --db="tmp.duckdb" \
            --parquet-dir="./parq" \
            --db-tmp-dir-path="./db_tmp/" \
            --no-multiallelics \
            --config="/app/config.json" \
            --batch-size=~{batch_size}
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

task merge_dbs {
    input {
        Array[Array[File]] duckdbs

        String docker_image = "us-central1-docker.pkg.dev/depmap-omics/terra-images/vcf-to-duckdb"
        String docker_image_hash_or_tag = ":production"
        Int mem_gb = 8
        Int cpu = 2
        Int preemptible = 1
        Int max_retries = 1
        Int additional_disk_gb = 0
    }

    Array[File] all_db_files = flatten(duckdbs)

    Int disk_space = ceil(3 * size(all_db_files, "GiB")) + 10 + additional_disk_gb

    command <<<
        set -euo pipefail

        # collect unique enclosing folder names of sharded DB files
        folder_names_file=$(mktemp)

        for f in ~{sep=' ' all_db_files}; do
          dirname "$f"
        done | sort -u > "$folder_names_file"

        # build repeated --db options
        db_opts=""

        while IFS= read -r path; do
            db_opts="$db_opts --db=$path"
        done < "$folder_names_file"

        python -m vcf_to_duckdb merge $db_opts --parquet-dir="./parq"
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
