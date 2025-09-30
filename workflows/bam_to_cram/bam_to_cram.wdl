version 1.0

workflow bam_to_cram {
    input {
        String sample_id
        File bam
        File bai
        File ref_fasta
        File ref_fasta_index
    }

    call convert_bam_to_cram {
        input:
            sample_id = sample_id,
            bam = bam,
            bai = bai,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index
    }

    output {
        File cram = convert_bam_to_cram.cram
        File crai = convert_bam_to_cram.crai
    }
}

task convert_bam_to_cram {
    input {
        String sample_id
        File bam
        File bai
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
        ceil(1.5 * size(bam, "GiB") + size(ref_fasta, "GiB")) + 10 + additional_disk_gb
    )

    Int n_threads = cpu

    command <<<
        set -euo pipefail

        echo "Converting BAM to CRAM"
        samtools view \
            -@ ~{n_threads} \
            -C \
            -T "~{ref_fasta}" \
            -o "~{sample_id}.cram" \
            "~{bam}"

        echo "Indexing CRAM"
        samtools index -@ ~{n_threads} "~{sample_id}.cram"

        mv "~{sample_id}.cram.crai" "~{sample_id}.crai"
    >>>

    output {
        File cram = "~{sample_id}.cram"
        File crai = "~{sample_id}.crai"
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
