version 1.0

workflow cram_to_bam {
    input {
        String sample_id
        File cram
        File crai
        File ref_fasta
        File ref_fasta_index
    }

    call convert_cram_to_bam {
        input:
            sample_id = sample_id,
            cram = cram,
            crai = crai,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index
    }

    output {
        File bam = convert_cram_to_bam.bam
        File bai = convert_cram_to_bam.bai
    }
}

task convert_cram_to_bam {
    input {
        String sample_id
        File cram
        File crai
        File ref_fasta
        File ref_fasta_index

        String docker_image = "us-central1-docker.pkg.dev/depmap-omics/terra-images/samtools"
        String docker_image_hash_or_tag = ":production"
        Int mem_gb = 4
        Int cpu = 16
        Int preemptible = 1
        Int max_retries = 0
        Int additional_disk_gb = 0
    }

    Int disk_space = (
        ceil(6 * size(cram, "GiB") + size(ref_fasta, "GiB")) + 10 + additional_disk_gb
    )

    Int n_threads = cpu

    command <<<
        set -euo pipefail

        echo "Converting CRAM to BAM"
        samtools view \
            -@ ~{n_threads} \
            -b \
            -T "~{ref_fasta}" \
            -o "~{sample_id}.bam" \
            "~{cram}"

        echo "Indexing BAM"
        samtools index -@ ~{n_threads} "~{sample_id}.bam"

        mv "~{sample_id}.bam.bai" "~{sample_id}.bai"
    >>>

    output {
        File bam = "~{sample_id}.bam"
        File bai = "~{sample_id}.bai"
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
