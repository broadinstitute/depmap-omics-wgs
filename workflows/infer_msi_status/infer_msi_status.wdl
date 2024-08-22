version 1.0

workflow infer_msi_status {
    input {
        String workflow_version = "1.0"
        String workflow_source_url # populated automatically with URL of this script

        String sample_id
        File bam
        File bai
    }

    call msisensor2 {
        input:
            sample_id=sample_id,
            bam=bam,
            bai=bai
    }

    output {
        Float msisensor2_score = msisensor2.msisensor2_score
        File msisensor2_output = msisensor2.msisensor2_output
        File msisensor2_output_dis = msisensor2.msisensor2_output_dis
        File msisensor2_output_somatic = msisensor2.msisensor2_output_somatic
    }
}

task msisensor2 {
    input {
        String sample_id
        File bam
        File bai

        Int preemptible = 2
        Int max_retries = 1
        Int cpu = 1
        Int mem_gb = 6
        Int additional_disk_gb = 0
    }

    String bam_path = basename(bam)
    Int disk_space = ceil(size(bam, "GiB")) + 10 + additional_disk_gb

    command <<<
        set -euo pipefail

        mv ~{bam} .
        mv ~{bai} .

        msisensor2 msi \
            -M /msisensor2/models_hg38 \
            -t ~{bam_path} \
            -o ~{sample_id}.msisensor2.output

        head -2 ~{sample_id}.msisensor2.output | \
            tail -1 | \
            cut -f3 > \
            ~{sample_id}.msisensor2.score
    >>>

    output {
        Float msisensor2_score = read_float("~{sample_id}.msisensor2.score")
        File msisensor2_output = "~{sample_id}.msisensor2.output"
        File msisensor2_output_dis = "~{sample_id}.msisensor2.output_dis"
        File msisensor2_output_somatic = "~{sample_id}.msisensor2.output_somatic"
    }

    runtime {
        docker: "us-docker.pkg.dev/depmap-omics/public/msisensor2:production"
        memory: "~{mem_gb} GiB"
        disks: "local-disk ~{disk_space} SSD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }
}
