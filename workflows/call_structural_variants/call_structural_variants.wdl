version 1.0

workflow call_structural_variants {
    input {
        String sample_id
        File tumor_bam
        File tumor_bai
        File? normal_bam
        File? normal_bai
        File ref_fasta
        File ref_fasta_index
        Boolean is_major_contigs_only = true
        File major_contig_bed
        File major_contig_bed_index
    }

    call run_manta {
        input: sample_id = sample_id,
            tumor_bam = tumor_bam,
            tumor_bai = tumor_bai,
            normal_bam = normal_bam,
            normal_bai = normal_bai,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            major_contig_bed = major_contig_bed,
            major_contig_bed_index = major_contig_bed_index,
            is_major_contigs_only = is_major_contigs_only
    }

    output {
        File? sv_diploid_vcf = run_manta.diploid_vcf
        File? sv_diploid_vcf_index = run_manta.diploid_vcf_index
        File sv_somatic_vcf = run_manta.somatic_vcf
        File sv_somatic_vcf_index = run_manta.somatic_vcf_index
        File sv_candidate_vcf = run_manta.candidate_vcf
        File sv_candidate_vcf_index = run_manta.candidate_vcf_index
        File sv_candidate_indel_vcf = run_manta.candidate_indel_vcf
        File sv_candidate_indel_vcf_index = run_manta.candidate_indel_vcf_index
    }
}

task run_manta {
    input {
        String sample_id
        File tumor_bam
        File tumor_bai
        File? normal_bam
        File? normal_bai
        File ref_fasta
        File ref_fasta_index
        Boolean is_major_contigs_only
        File major_contig_bed
        File major_contig_bed_index

        String docker_image = "us-central1-docker.pkg.dev/depmap-omics/terra-images/manta"
        String docker_image_hash_or_tag = ":production"
        Int preemptible = 1
        Int max_retries = 0
        Int cpu = 8
        Float mem_per_job_gb = 0.4
        Int additional_disk_gb = 0
    }

    Float jobs_per_cpu = 1.3
    Int num_jobs = round(cpu * jobs_per_cpu)
    Int mem_gb = ceil(num_jobs * mem_per_job_gb)
    Int disk_space = (
        ceil(size(tumor_bam, "GiB") + size(normal_bam, "GiB") + size(ref_fasta, "GiB"))
        + 10 + additional_disk_gb
     )

    command <<<
        set -euo pipefail

        touch "~{tumor_bai}"

        # link localized files to working dir
        ln -vs "~{tumor_bam}" "~{basename(tumor_bam)}"
        ln -vs "~{tumor_bai}" "~{basename(tumor_bai)}"

        if [[ -f "~{normal_bam}" ]]; then
            touch "~{normal_bai}"
            ln -vs "~{normal_bam}" "normal.bam"
            ln -vs "~{normal_bai}" "normal.bai"
            normal_command_line="--normalBam normal.bam"
        else
            normal_command_line=""
        fi

        ln -vs "~{ref_fasta}" "~{basename(ref_fasta)}"
        ln -vs "~{ref_fasta_index}" "~{basename(ref_fasta_index)}"

        if [ "~{is_major_contigs_only}" = "true" ]
        then
            ln -vs "~{major_contig_bed}" "~{basename(major_contig_bed)}"
            ln -vs "~{major_contig_bed_index}" "~{basename(major_contig_bed_index)}"
            major_contig_line="--callRegions ~{basename(major_contig_bed)}"
        else
            major_contig_line=""
        fi

        configManta.py \
            --tumorBam="~{basename(tumor_bam)}" \
            $normal_command_line \
            $major_contig_line \
            --referenceFasta="~{basename(ref_fasta)}" \
            --runDir="."
            --retainTempFiles

        # always tell manta there are 2 GiB per job, otherwise it will scale back the
        # requested number of jobs, even if they won't need that much memory
        ./runWorkflow.py \
            --mode="local" \
            --jobs=~{cpu} \
            --memGb=$((~{num_jobs} * 2)) \
            --quiet

        if [[ -f "~{normal_bam}" ]]; then
           mv results/variants/diploidSV.vcf.gz "~{sample_id}.diploid_sv.vcf.gz"
           mv results/variants/diploidSV.vcf.gz.tbi "~{sample_id}.diploid_sv.vcf.gz.tbi"
           mv results/variants/somaticSV.vcf.gz "~{sample_id}.somatic_sv.vcf.gz"
           mv results/variants/somaticSV.vcf.gz.tbi "~{sample_id}.somatic_sv.vcf.gz.tbi"
        else
           touch "~{sample_id}.diploid_sv.vcf.gz"
           touch "~{sample_id}.diploid_sv.vcf.gz.tbi"
           mv results/variants/tumorSV.vcf.gz "~{sample_id}.somatic_sv.vcf.gz"
           mv results/variants/tumorSV.vcf.gz.tbi "~{sample_id}.somatic_sv.vcf.gz.tbi"
        fi

        mv results/variants/candidateSV.vcf.gz "~{sample_id}.candidate_sv.vcf.gz"
        mv results/variants/candidateSV.vcf.gz.tbi "~{sample_id}.candidate_sv.vcf.gz.tbi"
        mv results/variants/candidateSmallIndels.vcf.gz "~{sample_id}.candidate_small_indels.vcf.gz"
        mv results/variants/candidateSmallIndels.vcf.gz.tbi "~{sample_id}.candidate_small_indels.vcf.gz.tbi"
    >>>

    output {
        File? diploid_vcf = "~{sample_id}.diploid_sv.vcf.gz"
        File? diploid_vcf_index = "~{sample_id}.diploid_sv.vcf.gz.tbi"
        File somatic_vcf = "~{sample_id}.somatic_sv.vcf.gz"
        File somatic_vcf_index = "~{sample_id}.somatic_sv.vcf.gz.tbi"
        File candidate_vcf = "~{sample_id}.candidate_sv.vcf.gz"
        File candidate_vcf_index = "~{sample_id}.candidate_sv.vcf.gz.tbi"
        File candidate_indel_vcf = "~{sample_id}.candidate_small_indels.vcf.gz"
        File candidate_indel_vcf_index = "~{sample_id}.candidate_small_indels.vcf.gz.tbi"
    }

    runtime {
        docker: "~{docker_image}~{docker_image_hash_or_tag}"
        memory: "~{mem_gb} GiB"
        disks: "local-disk ~{disk_space} SSD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }

    meta {
        allowNestedInputs: true
    }
}
