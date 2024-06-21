version 1.0

workflow annotate_mutations {
    input {
        String workflow_version = "1.0"
        String workflow_url # populate this with the public URL of this script

        File ref_fasta
        File ref_fasta_index
        File ref_dict
        String sample_id
        File input_vcf

        # scatter/gather VCF by chromosome
        File? xy_intervals

        # misc. bcftools-based annotation
        Boolean fix_ploidy = true
        Boolean filter_vcf = true
        Boolean annot_seg_dups = true
        Boolean annot_repeat_masker = true
        Boolean annot_hess_drivers = true
        Boolean normalize_indels = true
        String? exclude_string
        File? segdup_bed
        File? segdup_bed_index
        File? repeatmasker_bed
        File? repeatmasker_bed_index
        File? hess_drivers
        File? hess_drivers_index

        # snpEff and SnpSift annotation
        Boolean annot_snpeff = true
        Boolean annot_snpsift = true
        File? clinvar_vcf
        File? clinvar_vcf_index

        # Ensembl VEP annotation
        Boolean annot_ensembl_vep = true
        File? ref_fasta_bgz
        File? gene_constraint_scores
        File? loftool_scores
        File? alpha_missense
        File? alpha_missense_index
        String? vep_chrom_cache_url_prefix
        String? vep_chrom_cache_url_suffix

        # Funcotator annotation
        Boolean annot_funcotator = true
        String? reference_version = "hg38"
        Boolean? use_gnomad = true
        Boolean? filter_funcotations = false
        File? interval_list
        String? transcript_selection_mode
        Array[String]? transcript_selection_list
        Array[String]? funcotator_annotation_defaults
        Array[String]? funcotator_annotation_overrides
        String? funcotator_data_sources_url
        String? funcotator_extra_args
        String? gcs_project_for_requester_pays
    }

    call compress_vcf {
        input:
            vcf = input_vcf,
            output_file_base_name = sample_id
    }

    if (fix_ploidy || filter_vcf || annot_seg_dups || annot_repeat_masker || annot_hess_drivers || normalize_indels) {
        call annot_with_bcftools {
            input:
                vcf = compress_vcf.vcf_compressed,
                output_file_base_name = sample_id + "_bcftools_annot",
                fix_ploidy = fix_ploidy,
                filter_vcf = filter_vcf,
                annot_seg_dups = annot_seg_dups,
                annot_repeat_masker = annot_repeat_masker,
                annot_hess_drivers = annot_hess_drivers,
                normalize_indels = normalize_indels,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                exclude_string = exclude_string,
                segdup_bed = segdup_bed,
                segdup_bed_index = segdup_bed_index,
                repeatmasker_bed = repeatmasker_bed,
                repeatmasker_bed_index = repeatmasker_bed_index,
                hess_drivers = hess_drivers,
                hess_drivers_index = hess_drivers_index
        }
    }

    if (annot_snpeff || annot_snpsift) {
        call snpeff_snpsift {
            input:
                vcf = select_first([annot_with_bcftools.vcf_annot, compress_vcf.vcf_compressed]),
                annot_snpeff = annot_snpeff,
                annot_snpsift = annot_snpsift,
                output_file_base_name = sample_id + "_snpeff_snpsift_annot",
                clinvar_vcf = clinvar_vcf,
                clinvar_vcf_index = clinvar_vcf_index
        }
    }

    if (annot_ensembl_vep) {
        call split_vcf_by_chrom {
            input:
                vcf = select_first([
                    snpeff_snpsift.vcf_annot,
                    annot_with_bcftools.vcf_annot,
                    compress_vcf.vcf_compressed
                ]),
                xy_intervals = select_first([xy_intervals])
        }

        scatter (vcf in split_vcf_by_chrom.vcfs) {
            String chrom_num = sub(sub(basename(vcf), "^chr", ""), ".vcf.gz$", "")
            File vep_cache = vep_chrom_cache_url_prefix + chrom_num + vep_chrom_cache_url_suffix

            call ensembl_vep {
                input:
                    vcf = vcf,
                    output_file_base_name = sample_id + "_chr" + chrom_num + "_ensembl_vep_annot",
                    vep_cache = vep_cache,
                    ref_fasta_bgz = select_first([ref_fasta_bgz]),
                    gene_constraint_scores = select_first([gene_constraint_scores]),
                    loftool_scores = select_first([loftool_scores]),
                    alpha_missense = select_first([alpha_missense]),
                    alpha_missense_index = select_first([alpha_missense_index])
            }
        }

        call gather_vcfs {
            input:
                vcfs = ensembl_vep.vcf_annot,
                output_file_base_name = sample_id + "_ensembl_vep_annot",
        }
    }

    if (annot_funcotator) {
        call funcotate {
            input:
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                ref_dict = ref_dict,
                input_vcf = select_first([
                    gather_vcfs.output_vcf,
                    snpeff_snpsift.vcf_annot,
                    annot_with_bcftools.vcf_annot,
                    compress_vcf.vcf_compressed
                ]),
                reference_version = select_first([reference_version]),
                output_file_base_name = sample_id + "_funco_annot",
                output_format = "VCF",
                compress = true,
                use_gnomad = select_first([use_gnomad]),
                filter_funcotations = select_first([filter_funcotations]),
                funcotator_data_sources_url = select_first([funcotator_data_sources_url]),
                interval_list = select_first([interval_list]),
                transcript_selection_mode = transcript_selection_mode,
                transcript_selection_list = transcript_selection_list,
                funcotator_annotation_defaults = funcotator_annotation_defaults,
                funcotator_annotation_overrides = funcotator_annotation_overrides,
                extra_args = funcotator_extra_args,
                gcs_project_for_requester_pays = gcs_project_for_requester_pays
        }
    }

    File vcf_annot = select_first([
        funcotate.funcotated_output_file,
        gather_vcfs.output_vcf,
        snpeff_snpsift.vcf_annot,
        annot_with_bcftools.vcf_annot,
        compress_vcf.vcf_compressed
    ])

    call index_vcf {
        input:
            vcf = vcf_annot
    }

    output {
        File mut_annot_vcf = select_first([
            funcotate.funcotated_output_file,
            gather_vcfs.output_vcf,
            snpeff_snpsift.vcf_annot,
            annot_with_bcftools.vcf_annot,
            compress_vcf.vcf_compressed
        ])
        File mut_annot_vcf_index = index_vcf.vcf_index
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

task index_vcf {
    input {
        File vcf
        String index_format = "CSI"

        String docker_image
        String docker_image_hash_or_tag
        Int mem_gb = 4
        Int cpu = 1
        Int preemptible = 3
        Int max_retries = 0
        Int additional_disk_gb = 0
    }

    Int disk_space = ceil(size(vcf, "GiB")) + 10 + additional_disk_gb

    String vcf_basename = basename(vcf)
    String index_format_option = if index_format == "CSI" then "--csi" else "--tbi"
    String index_file_ext = if index_format == "CSI" then ".csi" else ".tbi"

    command <<<
        set -euo pipefail

        bcftools index \
            "~{vcf}" \
            --output="~{vcf_basename}~{index_file_ext}"
            ~{index_format_option} \
            --no-version
    >>>

    output {
        File vcf_index = "~{vcf_basename}~{index_file_ext}"
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

task split_vcf_by_chrom {
    input {
        File vcf
        File xy_intervals

        String docker_image
        String docker_image_hash_or_tag
        Int mem_gb = 2
        Int cpu = 1
        Int preemptible = 3
        Int max_retries = 0
        Int additional_disk_gb = 0
    }

    Int disk_space = ceil(2 * size(vcf, "GiB")) + 10 + additional_disk_gb

    command <<<
        set -euo pipefail

        bcftools index ~{vcf}

        mkdir vcfs

        for chr in $(cat ~{xy_intervals}); do
            split_out="vcfs/${chr}.vcf.gz"

            bcftools view \
                "~{vcf}" \
                --regions="${chr}" \
                --output="${split_out}" \
                --no-version
        done
    >>>

    output {
        Array[File] vcfs = glob("vcfs/*.vcf.gz")
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

task gather_vcfs {
    input {
        Array[File] vcfs
        String output_file_base_name

        String docker_image
        String docker_image_hash_or_tag
        Int mem_gb = 8
        Int cpu = 1
        Int preemptible = 3
        Int max_retries = 0
        Int additional_disk_gb = 0
    }

    Int command_mem_mb = 1000 * mem_gb - 500
    Int disk_space = ceil(3 * size(vcfs, "GiB")) + 10 + additional_disk_gb

    parameter_meta {
        vcfs: { localization_optional: true }
    }

    command <<<
        set -euo pipefail

        gatk --java-options "-Xmx~{command_mem_mb}m" \
        GatherVcfsCloud \
            --input ~{sep=" --input " vcfs} \
            --output "combined.vcf.gz" \
            --gather-type "BLOCK" \
            --disable-contig-ordering-check

        gatk --java-options "-Xmx~{command_mem_mb}m" \
            SortVcf \
            --INPUT "combined.vcf.gz" \
            --OUTPUT "~{output_file_base_name}.vcf.gz"
    >>>

    output {
        File output_vcf = "~{output_file_base_name}.vcf.gz"
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

task annot_with_bcftools {
    input {
        File vcf
        String output_file_base_name
        Boolean fix_ploidy
        Boolean filter_vcf
        Boolean annot_seg_dups
        Boolean annot_repeat_masker
        Boolean annot_hess_drivers
        Boolean normalize_indels
        File? ref_fasta
        File? ref_fasta_index
        String? exclude_string
        File? segdup_bed
        File? segdup_bed_index
        File? repeatmasker_bed
        File? repeatmasker_bed_index
        File? hess_drivers
        File? hess_drivers_index

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

        if ~{fix_ploidy}; then
            echo "Fixing ploidy"
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
        fi

        if ~{filter_vcf}; then
            echo "Filtering VCF: ~{exclude_string}"
            bcftools view \
                "~{vcf}" \
                --exclude='~{exclude_string}' \
                --output="~{vcf}.2"
            rm "~{vcf}" && mv "~{vcf}.2" "~{vcf}"
        fi

        if ~{annot_seg_dups}; then
            echo "Annotating segmental duplication regions"
            echo '##INFO=<ID=SEGDUP,Number=1,Type=String,Description="If variant is in a segmental duplication region">' > \
                segdup.hdr.vcf
            bcftools annotate \
                "~{vcf}" \
                --annotations="~{segdup_bed}" \
                --columns="CHROM,FROM,TO,SEGDUP" \
                --header-lines="segdup.hdr.vcf" \
                --output="~{vcf}.2"
            rm "~{vcf}" && mv "~{vcf}.2" "~{vcf}"
        fi

        if ~{annot_repeat_masker}; then
            echo "Annotating repeat masker regions"
            echo '##INFO=<ID=RM,Number=1,Type=String,Description="If variant is in a Repeat Masker region">' > \
                repeatmasker.hdr.vcf
            bcftools annotate \
                "~{vcf}" \
                --annotations="~{repeatmasker_bed}" \
                --columns="CHROM,FROM,TO,RM" \
                --header-lines="repeatmasker.hdr.vcf" \
                --output="~{vcf}.2"
            rm "~{vcf}" && mv "~{vcf}.2" "~{vcf}"
        fi

        if ~{annot_hess_drivers}; then
            echo "Annotating Hess drivers"
            echo '##INFO=<ID=HESS,Number=1,Type=String,Description="Hess driver signature">' > \
                hess.hdr.vcf
            bcftools annotate \
                "~{vcf}" \
                --annotations="~{hess_drivers}" \
                --columns="CHROM,POS,REF,ALT,HESS" \
                --header-lines="hess.hdr.vcf" \
                --output="~{vcf}.2"
            rm "~{vcf}" && mv "~{vcf}.2" "~{vcf}"
        fi

        if ~{normalize_indels}; then
            echo "Normalizing indels"
            bcftools norm "~{vcf}" \
                --multiallelics=- \
                --site-win=10000 \
                --fasta-ref="~{ref_fasta}" \
                --output="~{vcf}.2"
            rm "~{vcf}" && mv "~{vcf}.2" "~{vcf}"
        fi

        # ensure it's bgzipped
        bcftools view \
            "~{vcf}" \
            --output="~{output_file_base_name}.vcf.gz" \
            --no-version
    >>>

    output {
        File vcf_annot = "~{output_file_base_name}.vcf.gz"
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
        File vcf
        Boolean annot_snpeff
        Boolean annot_snpsift
        String output_file_base_name
        File? clinvar_vcf
        File? clinvar_vcf_index

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
        set -euo pipefail

        if ~{annot_snpeff}; then
            echo "Annotating with snpEff"
            java -Xmx14g -jar /app/snpEff.jar \
                ann \
                -noStats \
                GRCh38.mane.1.2.ensembl \
                "~{vcf}" > \
                "snpeff_out.vcf"

            bgzip "snpeff_out.vcf" -o "snpeff_out.vcf.gz"
            rm "snpeff_out.vcf" && mv "snpeff_out.vcf.gz" "~{vcf}"
        fi

        if ~{annot_snpsift}; then
            echo "Annotating with SnpSift"
            java -Xmx14g -jar /app/SnpSift.jar \
                annotate \
                -tabix \
                -noDownload \
                ~{clinvar_vcf} \
                "~{vcf}" > \
                "snpsift_out.vcf"

            bgzip "snpsift_out.vcf" -o "snpsift_out.vcf.gz"
            rm "snpsift_out.vcf" && mv "snpsift_out.vcf.gz" "~{vcf}"
        fi

        mv "~{vcf}" "~{output_file_base_name}.vcf.gz"
    >>>

    output {
        File vcf_annot = "~{output_file_base_name}.vcf.gz"
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

task ensembl_vep {
    input {
        File vcf
        String output_file_base_name
        File vep_cache
        File ref_fasta_bgz
        File gene_constraint_scores
        File loftool_scores
        File alpha_missense
        File alpha_missense_index

        String docker_image
        String docker_image_hash_or_tag
        Int mem_gb = 8
        Int cpu = 2
        Int preemptible = 3
        Int max_retries = 0
        Int additional_disk_gb = 0
    }

    Int disk_space = (
        ceil(10 * size(vcf, "GiB") +
             size(ref_fasta_bgz, "GiB") +
             size(gene_constraint_scores, "GiB") +
             size(loftool_scores, "GiB") +
             size(alpha_missense, "GiB") +
             3 * size(vep_cache, "GiB")
        ) + 10 + additional_disk_gb
    )

    command <<<
        set -euo pipefail

        gunzip -c "~{gene_constraint_scores}" | \
            awk '{ print $2, $20 }' > \
            "pLI.tsv"

        echo "Populating VEP cache data"
        mkdir -p "vep_cache"
        tar -C "vep_cache" -xzf "~{vep_cache}"

        echo "Running VEP"
        /opt/vep/src/ensembl-vep/vep \
            --species="homo_sapiens" \
            --assembly="GRCh38" \
            --dir_cache="vep_cache" \
            --dir_plugins="/plugins" \
            --input_file="~{vcf}" \
            --output_file="~{output_file_base_name}.vcf" \
            --plugin="pLI,pLI.tsv" \
            --plugin="LoFtool,~{loftool_scores}" \
            --plugin="AlphaMissense,file=~{alpha_missense}" \
            --fasta="~{ref_fasta_bgz}" \
            --fork=~{cpu} \
            --buffer_size=5000 \
            --offline \
            --cache \
            --vcf \
            --no_stats \
            --everything \
            --pick

        bgzip "~{output_file_base_name}.vcf" --threads=~{cpu}
    >>>

    output {
        File vcf_annot = "~{output_file_base_name}.vcf.gz"
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

# Adapted from https://github.com/broadinstitute/gatk/blob/master/scripts/funcotator_wdl/funcotator.wdl
#
# Modifications:
#   - Localize pre-extracted Funcotator datasource from GCS
#   - Use SSD instead of HDD by default
#
# Use this file as a base and manually implement changes in the upstream code in order
# to retain these modifications.
task funcotate {
    input {
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        String funcotator_data_sources_url

        File input_vcf

        String reference_version

        String output_file_base_name
        String output_format

        Boolean compress
        Boolean use_gnomad

        String? control_id
        String? case_id
        String? sequencing_center
        String? sequence_source
        String? transcript_selection_mode
        Array[String]? transcript_selection_list
        Array[String]? funcotator_annotation_defaults
        Array[String]? funcotator_annotation_overrides
        Array[String]? funcotator_excluded_fields
        Boolean? filter_funcotations
        File? interval_list

        String? extra_args
        String? gcs_project_for_requester_pays

        String docker_image
        String docker_image_hash_or_tag
        Int mem_gb = 3
        Int cpu = 1
        Int preemptible = 3
        Int max_retries = 0
        Int additional_disk_gb = 0
    }

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = mem_gb * 1000
    Int command_mem = machine_mem - 1000

    # Calculate disk size:
    Float ref_size_gb = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + size(ref_dict, "GiB")
    Float vcf_size_gb = size(input_vcf, "GiB")
    Float datasources_size_gb = 25
    Int disk_space = ceil(ref_size_gb + datasources_size_gb + vcf_size_gb) + 20

    # Process input args:
    String output_maf = output_file_base_name + ".maf"
    String output_maf_index = output_maf + ".idx"
    String output_vcf = output_file_base_name + if compress then ".vcf.gz" else ".vcf"
    String output_vcf_index = output_vcf +  if compress then ".tbi" else ".idx"
    String output_file = if output_format == "MAF" then output_maf else output_vcf
    String output_file_index = if output_format == "MAF" then output_maf_index else output_vcf_index
    String transcript_selection_arg = if defined(transcript_selection_list) then " --transcript-list " else ""
    String annotation_def_arg = if defined(funcotator_annotation_defaults) then " --annotation-default " else ""
    String annotation_over_arg = if defined(funcotator_annotation_overrides) then " --annotation-override " else ""
    String filter_funcotations_args = if defined(filter_funcotations) && (filter_funcotations) then " --remove-filtered-variants " else ""
    String excluded_fields_args = if defined(funcotator_excluded_fields) then " --exclude-field " else ""
    String interval_list_arg = if defined(interval_list) then " -L " else ""
    String extra_args_arg = select_first([extra_args, ""])

    parameter_meta {
        ref_fasta: { localization_optional: true }
        ref_fasta_index: { localization_optional: true }
        ref_dict: { localization_optional: true }
        input_vcf: { localization_optional: true }
    }

    command <<<
        set -e
        export GATK_LOCAL_JAR="/root/gatk.jar"

        # Hack to validate our WDL inputs:
        #
        # NOTE: This happens here so that we don't waste time copying down the data sources if there's an error.
        if [[ "~{output_format}" != "MAF" ]] && [[ "~{output_format}" != "VCF" ]] ; then
            echo "ERROR: Output format must be MAF or VCF."
        fi

        # download and extract Funcotator data sources
        DATA_SOURCES_FOLDER="./funcotator_data_sources"
        gcloud storage rsync --recursive ~{funcotator_data_sources_url} $DATA_SOURCES_FOLDER

        if ~{use_gnomad} ; then
            echo "Enabling gnomAD..."
            for potential_gnomad_gz in gnomAD_exome.tar.gz gnomAD_genome.tar.gz ; do
                if [[ -f $DATA_SOURCES_FOLDER/$potential_gnomad_gz ]] ; then
                    cd $DATA_SOURCES_FOLDER
                    tar -zvxf $potential_gnomad_gz
                    cd -
                else
                    echo "ERROR: Cannot find gnomAD folder: $potential_gnomad_gz" 1>&2
                    false
                fi
            done
        fi

        # Run Funcotator:
        gatk --java-options "-Xmx~{command_mem}m" Funcotator \
            --data-sources-path $DATA_SOURCES_FOLDER \
            --ref-version ~{reference_version} \
            --output-file-format ~{output_format} \
            -R ~{ref_fasta} \
            -V ~{input_vcf} \
            -O ~{output_file} \
            ~{interval_list_arg} ~{default="" interval_list} \
            --annotation-default normal_barcode:~{default="Unknown" control_id} \
            --annotation-default tumor_barcode:~{default="Unknown" case_id} \
            --annotation-default Center:~{default="Unknown" sequencing_center} \
            --annotation-default source:~{default="Unknown" sequence_source} \
            ~{"--transcript-selection-mode " + transcript_selection_mode} \
            ~{transcript_selection_arg}~{default="" sep=" --transcript-list " transcript_selection_list} \
            ~{annotation_def_arg}~{default="" sep=" --annotation-default " funcotator_annotation_defaults} \
            ~{annotation_over_arg}~{default="" sep=" --annotation-override " funcotator_annotation_overrides} \
            ~{excluded_fields_args}~{default="" sep=" --exclude-field " funcotator_excluded_fields} \
            ~{filter_funcotations_args} \
            ~{extra_args_arg} \
            ~{"--gcs-project-for-requester-pays " + gcs_project_for_requester_pays}

        # Make sure we have a placeholder index for MAF files so this workflow doesn't fail:
        if [[ "~{output_format}" == "MAF" ]] ; then
            touch ~{output_maf_index}
        fi
    >>>

    runtime {
        docker: "~{docker_image}~{docker_image_hash_or_tag}"
        memory: "~{mem_gb} GB"
        disks: "local-disk ~{disk_space} SSD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }

    output {
        File funcotated_output_file = "~{output_file}"
    }
}
