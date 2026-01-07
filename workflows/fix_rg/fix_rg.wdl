version 1.0


workflow fix_rg {
    input {
        String sample_id
        File input_bam
    }

    call do_fix_rg {
        input:
            sample_id = sample_id,
            input_bam = input_bam
    }

    output {
        File bam_fixed = do_fix_rg.bam_fixed
        File bai_fixed = do_fix_rg.bai_fixed
    }
}

task do_fix_rg {
    input {
        String sample_id
        File input_bam

        String docker_image
        String docker_image_hash_or_tag
        Int mem_gb = 8
        Int cpu = 2
        Int preemptible = 0
        Int max_retries = 0
        Int additional_disk_gb = 0
    }

    Int disk_space = ceil(3 * size(input_bam, "GiB")) + additional_disk_gb

    command <<<
        samtools view -H ~{input_bam} \
            | awk 'BEGIN{OFS="\t"}
                   /^@RG/ {
                     sm=""
                     for (i=1;i<=NF;i++) {
                       if ($i ~ /^SM:/) sm=$i
                     }
                     if (sm != "") print $0, "LB:" substr(sm,4)
                     else print
                     next
                   }
                   {print}' \
            > header.sam

        samtools reheader header.sam ~{input_bam} > ~{sample_id}.fixed.bam
        samtools index ~{sample_id}.fixed.bam -o ~{sample_id}.fixed.bai
    >>>

    output {
        File bam_fixed = "~{sample_id}.fixed.bam"
        File bai_fixed = "~{sample_id}.fixed.bai"
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
