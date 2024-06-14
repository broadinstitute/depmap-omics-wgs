version 1.0

workflow annotate_mutations {
    input {
        String workflow_version = "1.0"
        String workflow_url # populate this with the public URL of this script

        File ref_fasta
        File ref_fasta_index
        String sample_id
        File input_vcf
        File input_vcf_idx
        String exclude_string
        File segdup_bed
        File segdup_bed_index
        File repeatmasker_bed
        File repeatmasker_bed_index
        File clinvar_vcf
        File clinvar_vcf_index
        String output_file_base_name
        String reference_version = "hg38"
        String output_format = "VCF"
        Boolean compress = true
        Boolean use_gnomad = true

        File? interval_list
        String? transcript_selection_mode
        Array[String]? transcript_selection_list
        Array[String]? annotation_defaults
        Array[String]? annotation_overrides
        String funcotator_data_sources_url
        String? funcotator_extra_args
        String? gcs_project_for_requester_pays

        File? gatk_override
    }

    call fix_ploidy {
        input:
            sample_id = sample_id,
            vcf = input_vcf,
            vcf_idx = input_vcf_idx
    }

    call filter_and_mask {
        input:
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            sample_id = sample_id,
            vcf = fix_ploidy.vcf_fixedploidy,
            exclude_string = exclude_string,
            segdup_bed = segdup_bed,
            segdup_bed_index = segdup_bed_index,
            repeatmasker_bed = repeatmasker_bed,
            repeatmasker_bed_index = repeatmasker_bed_index
    }

    call snpeff_snpsift {
        input:
            sample_id = sample_id,
            vcf = filter_and_mask.vcf_filtered,
            clinvar_vcf = clinvar_vcf,
            clinvar_vcf_index = clinvar_vcf_index
    }

    output {
    }
}

task fix_ploidy {
    input {
        String sample_id
        File vcf
        File vcf_idx

        String docker_image
        String docker_image_hash_or_tag
        Int mem_gb = 4
        Int cpu = 1
        Int preemptible = 3
        Int max_retries = 0
        Int additional_disk_gb = 0
    }

    # we need twice the ungzipped size of the input VCF, plus the gzipped size of the
    # output VCF
    Int disk_space = ceil(2 * 10 * size(vcf, "GiB")) + ceil(size(vcf, "GiB")) + 10 + additional_disk_gb

    command <<<
        set -euo pipefail

        # set the genotype annotation to be homozygous when we have no ref reads and at
        # least 3 alt reads, e.g., 0/1/2 --> 1/2, 0/1/2/3 --> 1/2/3, etc.
        bcftools +setGT "~{vcf}" -- \
            -t q -i'INFO/DP>8 & AF>0.9' -n c:'m|m' > \
            "~{vcf}.2"
        bcftools +setGT "~{vcf}.2" -- \
            -t q -i'AD[*:0]=0 & INFO/DP>8 & GT="0/1/2"' -n c:'1/2' > \
            "~{vcf}"
        bcftools +setGT "~{vcf}" -- \
            -t q -i'AD[*:0]=0 & INFO/DP>8 & GT="0/1/2/3"' -n c:'1/2/3' > \
            "~{vcf}.2"
        bcftools +setGT "~{vcf}.2" -- \
            -t q -i'AD[*:0]=0 & INFO/DP>8 & GT="0/1/2/3/4"' -n c:'1/2/3/4' > \
            "~{vcf}"
        bcftools +setGT "~{vcf}" -- \
            -t q -i'AD[*:0]=0 & INFO/DP>8 & GT="0/1/2/3/4/5"' -n c:'1/2/3/4/5' > \
            "~{vcf}.2"
        bcftools +setGT "~{vcf}.2" -- \
            -t q -i'AD[*:0]=0 & INFO/DP>8 & GT="0/1/2/3/4/5/6"' -n c:'1/2/3/4/5/6' > \
            "~{vcf}"

        bcftools view "~{vcf}" -o "~{sample_id}_fixedploidy.vcf.gz"
    >>>

    output {
        File vcf_fixedploidy = "~{sample_id}_fixedploidy.vcf.gz"
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

task filter_and_mask {
    input {
        File ref_fasta
        File ref_fasta_index
        String sample_id
        File vcf
        String exclude_string
        File segdup_bed
        File segdup_bed_index
        File repeatmasker_bed
        File repeatmasker_bed_index

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

        echo "Filtering VCF: ~{exclude_string}"
        bcftools view \
            "~{vcf}" \
            --exclude='~{exclude_string}' \
            --output="~{sample_id}_filtered.vcf.gz" && rm "~{vcf}"

        echo "Annotating segmental duplication regions"
        echo '##INFO=<ID=SEGDUP,Number=1,Type=String,Description="If variant is in a segmental duplication region">' > \
            segdup.hdr.vcf
        bcftools annotate \
            "~{sample_id}_filtered.vcf.gz" \
            --annotations="~{segdup_bed}" \
            --columns="CHROM,FROM,TO,SEGDUP" \
            --header-lines="segdup.hdr.vcf" \
            --output="~{sample_id}_segdup.vcf.gz" && rm "~{sample_id}_filtered.vcf.gz"

        echo "Annotating repeat masker regions"
        echo '##INFO=<ID=RM,Number=1,Type=String,Description="If variant is in a Repeat Masker region">' > \
            repeatmasker.hdr.vcf
        bcftools annotate \
            "~{sample_id}_segdup.vcf.gz" \
            --annotations="~{repeatmasker_bed}" \
            --columns="CHROM,FROM,TO,RM" \
            --header-lines="repeatmasker.hdr.vcf" \
            --output="~{sample_id}_masked.vcf.gz" && rm "~{sample_id}_segdup.vcf.gz"

        echo "Normalizing indels"
        bcftools norm "~{sample_id}_masked.vcf.gz" \
            --multiallelics=- \
            --site-win=10000 \
            --fasta-ref="~{ref_fasta}" \
            --output="~{sample_id}_normalized.vcf.gz" && rm "~{sample_id}_masked.vcf.gz"

        echo "Reindexing VCF"
        bcftools index \
            "~{sample_id}_normalized.vcf.gz" \
            --csi \
            --output="~{sample_id}_normalized.vcf.gz.csi"
    >>>

    output {
        File vcf_filtered = "~{sample_id}_normalized.vcf.gz"
        File vcf_filtered_index = "~{sample_id}_normalized.vcf.gz.csi"
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

task snpeff_snpsift {
    input {
        String sample_id
        File vcf
        File clinvar_vcf
        File clinvar_vcf_index

        String docker_image
        String docker_image_hash_or_tag
        Int mem_gb = 16
        Int cpu = 1
        Int preemptible = 3
        Int max_retries = 0
        Int additional_disk_gb = 0
    }

    Int disk_space = ceil(2 * 10 * size(vcf, "GiB")) + 10 + additional_disk_gb

    command <<<
        echo "Annotating with snpEff"
        java -Xmx14g -jar /app/snpEff.jar \
            ann \
            -noStats \
            GRCh38.mane.1.2.ensembl \
            "~{vcf}" > \
            "~{sample_id}_snpeff.vcf" && rm "~{vcf}"

        echo "Annotating with SnpSift"
        java -Xmx14g -jar /app/SnpSift.jar \
            annotate \
            -tabix \
            -noDownload \
            ~{clinvar_vcf} \
            "~{sample_id}_snpeff.vcf" > \
            "~{sample_id}_snpsift.vcf" && rm "~{sample_id}_snpeff.vcf"

        bgzip "~{sample_id}_snpsift.vcf" -o "~{sample_id}_snpsift.vcf.gz"
    >>>

    output {
        File vcf_annot = "~{sample_id}_snpsift.vcf.gz"
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
