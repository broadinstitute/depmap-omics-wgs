version 1.0

workflow make_guide_mutation_beds {
    input {
        String sample_id
        File vcf
        File? vcf_idx
        File bed_chromium_v1
        File bed_chromium_v2
    }

    call filter_vcfs_by_bed {
        input:
            sample_id = sample_id,
            vcf = vcf,
            vcf_idx = vcf_idx,
            bed_chromium_v1 = bed_chromium_v1,
            bed_chromium_v2 = bed_chromium_v2
    }

    call intersect_beds {
        input:
            sample_id = sample_id,
            mut_bed_chromium_v1 = filter_vcfs_by_bed.mut_bed_chromium_v1,
            mut_bed_chromium_v2 = filter_vcfs_by_bed.mut_bed_chromium_v2,
            bed_chromium_v1 = bed_chromium_v1,
            bed_chromium_v2 = bed_chromium_v2
    }

    output {
        File guide_bed_chromium_v1 = intersect_beds.guide_bed_chromium_v1
        File guide_bed_chromium_v2 = intersect_beds.guide_bed_chromium_v2
    }
}

task filter_vcfs_by_bed {
    input {
        String sample_id
        File vcf
        File? vcf_idx
        File bed_chromium_v1
        File bed_chromium_v2

        String docker_image = "us-central1-docker.pkg.dev/depmap-omics/terra-images/bcftools"
        String docker_image_hash_or_tag = ":production"
        Int mem_gb = 4
        Int cpu = 1
        Int preemptible = 2
        Int max_retries = 0
        Int additional_disk_gb = 0
    }

    Int disk_space = ceil(3 * size(vcf, "GiB")) + 10 + additional_disk_gb

    command <<<
        set -euo pipefail

        if [[ ! -f "~{vcf_idx}" ]]; then
            bcftools index "~{vcf}"
        fi

        # create parallel arrays of labels and library BED files
        labels=("chromium_v1" "chromium_v2")
        beds=(~{bed_chromium_v1} ~{bed_chromium_v2})

        # loop over the tuples
        for i in "${!labels[@]}"; do
            label="${labels[$i]}"
            bed="${beds[$i]}"

            echo "Filtering VCF with ${bed}"

            bcftools query \
                --exclude="FILTER!='PASS'&GT!='mis'&GT!~'\.'" \
                --regions-file="$bed" \
                --format="%CHROM\t%POS\t%END\t%ALT{0}\n" \
                "~{vcf}" > "~{sample_id}.mut_${label}.bed"
        done
    >>>

    output {
        File mut_bed_chromium_v1 = "~{sample_id}.mut_chromium_v1.bed"
        File mut_bed_chromium_v2 = "~{sample_id}.mut_chromium_v2.bed"
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

task intersect_beds {
    input {
        String sample_id
        File mut_bed_chromium_v1
        File mut_bed_chromium_v2
        File bed_chromium_v1
        File bed_chromium_v2

        String docker_image = "us-central1-docker.pkg.dev/depmap-omics/terra-images/bedtools2"
        String docker_image_hash_or_tag = ":production"
        Int mem_gb = 4
        Int cpu = 1
        Int preemptible = 2
        Int max_retries = 0
        Int additional_disk_gb = 0
    }

    Int disk_space = ceil(
        size(mut_bed_chromium_v1, "GiB")
        + size(mut_bed_chromium_v2, "GiB")
    ) + 10 + additional_disk_gb

    command <<<
        set -euo pipefail

        # create parallel arrays of labels, mutation BED files, and library BED files
        labels=("chromium_v1" "chromium_v2")
        mut_beds=(~{mut_bed_chromium_v1} ~{mut_bed_chromium_v2})
        beds=(~{bed_chromium_v1} ~{bed_chromium_v2})

        # loop over the tuples
        for i in "${!labels[@]}"; do
            label="${labels[$i]}"
            mut_bed="${mut_beds[$i]}"
            bed="${beds[$i]}"

            echo "Intersecting ${mut_bed} with ${bed}"

            bedtools intersect \
                -a "${bed}" \
                -b "${mut_bed}" \
                -c \
            | bedtools sort -i - \
            > "~{sample_id}.guide_${label}.bed"
        done
    >>>

    output {
        File guide_bed_chromium_v1 = "~{sample_id}.guide_chromium_v1.bed"
        File guide_bed_chromium_v2 = "~{sample_id}.guide_chromium_v2.bed"
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
