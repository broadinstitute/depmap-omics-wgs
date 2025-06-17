version 1.0

workflow prep_annotations {
    input {
        String workflow_version = "1.0"
        String workflow_source_url = "" # populated automatically with URL of this script

        File ref_fasta
        File ref_fasta_index
        File ref_dict
        String sample_id
        File input_vcf

        # bcftools-based fixes and filters
        Boolean fix_ploidy = true
        Boolean filter_vcf = true
        Boolean normalize_indels = true
        String? exclude_string
    }

    call compress_vcf {
        input:
            vcf = input_vcf,
            output_file_base_name = sample_id
    }

    call fix_with_bcftools {
        input:
            vcf = compress_vcf.vcf_compressed,
            output_file_base_name = sample_id + "_bcftools_fixed",
            fix_ploidy = fix_ploidy,
            filter_vcf = filter_vcf,
            exclude_string = exclude_string,
            normalize_indels = normalize_indels,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index
    }

    output {
        File mut_annot_bcftools_fixed_vcf = fix_with_bcftools.vcf_fixed
    }
}

task compress_vcf {
    input {
        File vcf
        String output_file_base_name

        String docker_image
        String docker_image_hash_or_tag
        Int mem_gb = 4
        Int cpu = 1
        Int preemptible = 3
        Int max_retries = 0
        Int additional_disk_gb = 0
    }

    Int disk_space = ceil(2 * size(vcf, "GiB")) + 10 + additional_disk_gb

    command <<<
        set -euo pipefail

        bcftools view \
            "~{vcf}" \
            --output="~{output_file_base_name}.vcf.gz" \
            --no-version
    >>>

    output {
        File vcf_compressed = "~{output_file_base_name}.vcf.gz"
    }

    runtime {
        docker: "~{docker_image}~{docker_image_hash_or_tag}"
        memory: "~{mem_gb} GB"
        disks: "local-disk ~{disk_space} SSD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
        continueOnReturnCode: [0, 9] # FilterMutectCalls generates weird vcf.gz files
    }

    meta {
        allowNestedInputs: true
    }
}

task fix_with_bcftools {
    input {
        File vcf
        String output_file_base_name
        Boolean fix_ploidy
        Boolean filter_vcf
        Boolean normalize_indels
        File? ref_fasta
        File? ref_fasta_index
        String? exclude_string

        String docker_image
        String docker_image_hash_or_tag
        Int mem_gb = 4
        Int cpu = 1
        Int preemptible = 3
        Int max_retries = 0
        Int additional_disk_gb = 0
    }

    Int disk_space = (
        ceil(size(vcf, "GiB") + 2 * 10 * size(vcf, "GiB")) + 10 + additional_disk_gb
    )

    command <<<
        set -euo pipefail

        TMP_VCF="~{basename(vcf, '.vcf.gz')}_tmp.vcf.gz"

        # fix for https://github.com/broadinstitute/gatk/issues/6857
        # temporarily convert to text VCF
        echo "Fixing AS_FilterStatus allele delimiter"
        bcftools view "~{vcf}" -o tmp.vcf && rm "~{vcf}"

        awk '{
            # find the filter status field
            match($0, /\tAS_FilterStatus=[^;]+;/);

            if (RSTART != 0) {
                # the line matches
                before = substr($0, 1, RSTART-1);
                match_str = substr($0, RSTART, RLENGTH);
                after = substr($0, RSTART+RLENGTH);

                # temp replace "|" with a non-printing char to swap "|" and "," chars
                gsub(/\|/, "\x1e", match_str);
                gsub(/,/, "|", match_str);
                gsub(/\x1e/, ",", match_str);

                # print modified line
                print before match_str after;
            } else {
                # no match
                print $0;
            }
        }' tmp.vcf > swapped.vcf && rm tmp.vcf

        bcftools view swapped.vcf -o "~{vcf}" && rm swapped.vcf

        if ~{fix_ploidy}; then
            echo "Fixing ploidy"
            # set the genotype annotation to be homozygous when we have no ref reads and at
            # least 3 alt reads, e.g., 0/1/2 --> 1/2, 0/1/2/3 --> 1/2/3, etc.
            bcftools +setGT "~{vcf}" \
                -- -t q -i'INFO/DP>8 & AF>0.9' -n c:'m|m' \
                > "${TMP_VCF}"
            bcftools +setGT "${TMP_VCF}" \
                -- -t q -i'AD[*:0]=0 & INFO/DP>8 & GT="0/1/2"' -n c:'1/2' \
                > "~{vcf}"
            bcftools +setGT "~{vcf}" \
                -- -t q -i'AD[*:0]=0 & INFO/DP>8 & GT="0/1/2/3"' -n c:'1/2/3' \
                > "${TMP_VCF}"
            bcftools +setGT "${TMP_VCF}" \
                -- -t q -i'AD[*:0]=0 & INFO/DP>8 & GT="0/1/2/3/4"' -n c:'1/2/3/4' \
                > "~{vcf}"
            bcftools +setGT "~{vcf}" \
                -- -t q -i'AD[*:0]=0 & INFO/DP>8 & GT="0/1/2/3/4/5"' -n c:'1/2/3/4/5' \
                > "${TMP_VCF}"
            bcftools +setGT "${TMP_VCF}" \
                -- -t q -i'AD[*:0]=0 & INFO/DP>8 & GT="0/1/2/3/4/5/6"' -n c:'1/2/3/4/5/6' \
                > "~{vcf}"
        fi

        if ~{filter_vcf}; then
            echo "Filtering VCF: ~{exclude_string}"
            bcftools view \
                "~{vcf}" \
                --exclude='~{exclude_string}' \
                --output="${TMP_VCF}"
            rm "~{vcf}" && mv "${TMP_VCF}" "~{vcf}"
        fi

        if ~{normalize_indels}; then
            echo "Normalizing indels"
            bcftools norm "~{vcf}" \
                --multiallelics=- \
                --site-win=10000 \
                --fasta-ref="~{ref_fasta}" \
                --output="${TMP_VCF}"
            rm "~{vcf}" && mv "${TMP_VCF}" "~{vcf}"
        fi

        # ensure it's bgzipped
        bcftools view \
            "~{vcf}" \
            --output="~{output_file_base_name}.vcf.gz" \
            --no-version
    >>>

    output {
        File vcf_fixed = "~{output_file_base_name}.vcf.gz"
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
