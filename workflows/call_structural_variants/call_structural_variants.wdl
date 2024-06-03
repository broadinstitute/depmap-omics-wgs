version 1.0

workflow call_structural_variants {
    input {
        String workflow_version = "1.0"
        String workflow_url # populate this with the public URL of this script

        String sample_id
        File bam
        File bai
        File ref_fasta
        File ref_fasta_index
    }

    call run_manta {
        input: sample_id = sample_id,
            bam = bam,
            bai = bai,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index
    }

    output {
        File somatic_sv_vcf = run_manta.somatic_sv_vcf
        File somatic_sv_vcf_index = run_manta.somatic_sv_vcf_index
    }
}

task run_manta {
    input {
        String sample_id
        File bam
        File bai
        File ref_fasta
        File ref_fasta_index

        Int preemptible = 2
        Int max_retries = 1
        Int cpu = 8
        Int mem_per_job_gb = 0.4
        Int additional_disk_gb = 0
    }

    Float jobs_per_cpu = 1.3
    Int num_jobs = round(cpu * jobs_per_cpu)
    Float mem_gb = num_jobs * mem_per_job_gb
    Int disk_space = ceil(size(bam, "GiB")) + 10 + additional_disk_gb

    command {
        # link localized files to working dir
        ln -vs ~{bam} tumor.bam
        ln -vs ~{bai} tumor.bai
        ln -vs ~{ref_fasta} reference.fasta
        ln -vs ~{ref_fasta_index} reference.fasta.fai

        configManta.py \
            --tumorBam="tumor.bam" \
            --referenceFasta="reference.fasta" \
            --runDir="."

        # always tell manta there are 2 GiB per job, otherwise it will scale back the
        # requested number of jobs, even if they won't need that much memory
        ./runWorkflow.py \
            --mode="local" \
            --jobs=~{cpu}
            --memGb=$((~{num_jobs} * 2))

        mv results/variants/tumorSV.vcf.gz ~{sample_id}.somaticSV.vcf.gz
        mv results/variants/tumorSV.vcf.gz.tbi ~{sample_id}.somaticSV.vcf.gz.tbi
    }

    output {
        File somatic_sv_vcf = "${sample_id}.somaticSV.vcf.gz"
        File somatic_sv_vcf_index = "${sample_id}.somaticSV.vcf.gz.tbi"
    }

    runtime {
        docker: "us-docker.pkg.dev/depmap-omics/public/manta:production"
        memory: "~{mem_gb} GiB"
        disks: "local-disk ~{disk_space} SSD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }
}
