version 1.0

workflow make_guide_mutation_beds {
    input {
        String sample_id
        File vcf
        File? vcf_idx
        File bed_avana
        File bed_brunello
        File bed_humagne
        File bed_ky
        File bed_tko
        Float min_af
    }

    call filter_vcfs_by_bed {
        input:
            sample_id = sample_id,
            vcf = vcf,
            vcf_idx = vcf_idx,
            bed_avana = bed_avana,
            bed_brunello = bed_brunello,
            bed_humagne = bed_humagne,
            bed_ky = bed_ky,
            bed_tko = bed_tko,
            min_af = min_af
    }

    call intersect_beds {
        input:
            sample_id = sample_id,
            mut_bed_avana = filter_vcfs_by_bed.mut_bed_avana,
            mut_bed_brunello = filter_vcfs_by_bed.mut_bed_brunello,
            mut_bed_humagne = filter_vcfs_by_bed.mut_bed_humagne,
            mut_bed_ky = filter_vcfs_by_bed.mut_bed_ky,
            mut_bed_tko = filter_vcfs_by_bed.mut_bed_tko,
            bed_avana = bed_avana,
            bed_brunello = bed_brunello,
            bed_humagne = bed_humagne,
            bed_ky = bed_ky,
            bed_tko = bed_tko
    }

    output {
        File guide_bed_avana = intersect_beds.guide_bed_avana
        File guide_bed_brunello = intersect_beds.guide_bed_brunello
        File guide_bed_humagne = intersect_beds.guide_bed_humagne
        File guide_bed_ky = intersect_beds.guide_bed_ky
        File guide_bed_tko = intersect_beds.guide_bed_tko
    }
}

task filter_vcfs_by_bed {
    input {
        String sample_id
        File vcf
        File? vcf_idx
        File bed_avana
        File bed_brunello
        File bed_humagne
        File bed_ky
        File bed_tko
        Float min_af

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

        # ensure VCF is block gzipped, not plain gzipped (required for indexing)
        vcf_converted=false

        if [[ "~{vcf}" == *.gz ]]; then
            if ! file "~{vcf}" | grep -q "BGZF"; then
                echo "VCF is not block gzipped, converting..."
                bcftools view -O z -o temp.vcf.gz "~{vcf}"
                mv temp.vcf.gz "~{vcf}"
                vcf_converted=true
            fi
        fi

        if [[ ! -f "~{vcf_idx}" ]] || [[ "$vcf_converted" == "true" ]]; then
            # remove stale index if VCF was converted
            if [[ "$vcf_converted" == "true" ]] && [[ -f "~{vcf_idx}" ]]; then
                rm "~{vcf_idx}"
            fi

            # index the VCF now
            echo "Indexing VCF..."
            bcftools index "~{vcf}"
        fi

        # create parallel arrays of labels and library BED files
        labels=("avana" "brunello" "humagne" "ky" "tko")
        beds=(~{bed_avana} ~{bed_brunello} ~{bed_humagne} ~{bed_ky} ~{bed_tko})

        # loop over the tuples
        for i in "${!labels[@]}"; do
            label="${labels[$i]}"
            bed="${beds[$i]}"

            echo "Filtering VCF with ${bed}"

            bcftools query \
                --include="FILTER='PASS' & FORMAT/AF>=~{min_af}" \
                --regions-file="$bed" \
                --format="%CHROM\t%POS\t%END\t%ALT{0}[\t%AF]\n" \
                "~{vcf}" > "~{sample_id}.mut_${label}.bed"
        done
    >>>

    output {
        File mut_bed_avana = "~{sample_id}.mut_avana.bed"
        File mut_bed_brunello = "~{sample_id}.mut_brunello.bed"
        File mut_bed_humagne = "~{sample_id}.mut_humagne.bed"
        File mut_bed_ky = "~{sample_id}.mut_ky.bed"
        File mut_bed_tko = "~{sample_id}.mut_tko.bed"
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
        File mut_bed_avana
        File mut_bed_brunello
        File mut_bed_humagne
        File mut_bed_ky
        File mut_bed_tko
        File bed_avana
        File bed_brunello
        File bed_humagne
        File bed_ky
        File bed_tko

        String docker_image = "us-central1-docker.pkg.dev/depmap-omics/terra-images/bedtools2"
        String docker_image_hash_or_tag = ":production"
        Int mem_gb = 4
        Int cpu = 1
        Int preemptible = 2
        Int max_retries = 0
        Int additional_disk_gb = 0
    }

    Int disk_space = ceil(
        size(mut_bed_avana, "GiB")
        + size(mut_bed_brunello, "GiB")
        + size(mut_bed_humagne, "GiB")
        + size(mut_bed_ky, "GiB")
        + size(mut_bed_tko, "GiB")
    ) + 10 + additional_disk_gb

    command <<<
        set -euo pipefail

        # create parallel arrays of labels, mutation BED files, and library BED files
        labels=("avana" "brunello" "humagne" "ky" "tko")
        mut_beds=(~{mut_bed_avana} ~{mut_bed_brunello} ~{mut_bed_humagne} ~{mut_bed_ky} ~{mut_bed_tko})
        beds=(~{bed_avana} ~{bed_brunello} ~{bed_humagne} ~{bed_ky} ~{bed_tko})

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
        File guide_bed_avana = "~{sample_id}.guide_avana.bed"
        File guide_bed_brunello = "~{sample_id}.guide_brunello.bed"
        File guide_bed_humagne = "~{sample_id}.guide_humagne.bed"
        File guide_bed_ky = "~{sample_id}.guide_ky.bed"
        File guide_bed_tko = "~{sample_id}.guide_tko.bed"
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
