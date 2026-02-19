version 1.0

workflow compute_molecular_signatures {
    meta {
        description: "TODO"
    }

    parameter_meta {
        # inputs
        sample_id: "ID of this sample"

        # outputs
        todo: "TODO"
    }

    input {
        # per-sample inputs
        String sample_id

        # references
        String data_type
    }

    call step1 {
        input:
            sample_id = sample_id
    }

    call step1 {
        input:
            sample_id = sample_id
    }

    output {
    }
}

task step1 {
    meta {
        description: "TODO"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        sample_id: "ID of this sample"

        # outputs
        todo: "TODO"
    }

    input {
        String sample_id

        String docker_image  = "us-central1-docker.pkg.dev/depmap-omics/terra-images/compute-molecular-signatures"
        String docker_image_hash_or_tag = ":production"
        Int cpu = 2
        Int mem_gb = 16
        Int preemptible = 1
        Int max_retries = 0
        Int additional_disk_gb = 0
    }

    Int disk_space = 20 + additional_disk_gb

    command <<<
        set -euo pipefail

        python -m compute_molecular_signatures \
            filter-sample-gtf \
            --sample-id="~{sample_id}"
    >>>

    output {
        File todo = "~{sample_id}"
    }

    runtime {
        docker: "~{docker_image}~{docker_image_hash_or_tag}"
        memory: mem_gb + " GiB"
        disks: "local-disk " + disk_space + " SSD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }
}

task step2 {
    meta {
        description: "TODO"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        sample_id: "ID of this sample"

        # outputs
        todo: "TODO"
    }

    input {
        String sample_id

        String docker_image  = "us-central1-docker.pkg.dev/depmap-omics/terra-images/compute-molecular-signatures"
        String docker_image_hash_or_tag = ":production"
        Int cpu = 2
        Int mem_gb = 16
        Int preemptible = 1
        Int max_retries = 0
        Int additional_disk_gb = 0
    }

    Int disk_space = 20 + additional_disk_gb

    command <<<
        set -euo pipefail

        python -m compute_molecular_signatures \
            filter-sample-gtf \
            --sample-id="~{sample_id}"
    >>>

    output {
        File todo = "~{sample_id}"
    }

    runtime {
        docker: "~{docker_image}~{docker_image_hash_or_tag}"
        memory: mem_gb + " GiB"
        disks: "local-disk " + disk_space + " SSD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }
}
