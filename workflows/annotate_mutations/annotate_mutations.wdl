version 1.0

workflow annotate_mutations {
    input {
        String workflow_version = "1.0"
        String workflow_url # populate this with the public URL of this script

        String gatk_docker
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        File input_vcf
        File input_vcf_idx
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
        String? funcotator_extra_args
        String? gcs_project_for_requester_pays

        File? gatk_override
    }

    call Funcotate {
        input:
            gatk_docker = gatk_docker,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            input_vcf = input_vcf,
            input_vcf_idx = input_vcf_idx,
            reference_version = reference_version,
            output_file_base_name = output_file_base_name,
            output_format = output_format,
            compress = compress,
            use_gnomad = use_gnomad,

            interval_list = interval_list,
            transcript_selection_mode = transcript_selection_mode,
            transcript_selection_list = transcript_selection_list,
            annotation_defaults = annotation_defaults,
            annotation_overrides = annotation_overrides,
            extra_args = funcotator_extra_args,
            gcs_project_for_requester_pays = gcs_project_for_requester_pays,

            gatk_override = gatk_override
    }

    output {
        File funcotated_file_out = Funcotate.funcotated_output_file
        File funcotated_file_out_idx = Funcotate.funcotated_output_file_index
    }
}

# Adapted from https://github.com/broadinstitute/gatk/blob/master/scripts/funcotator_wdl/funcotator.wdl
#
# Modifications:
#   - Get latest Funcotator datasource archive dynamically
#   - Use SSD instead of HDD by default
#
# Use this file as a base and manually implement changes in the upstream code in order
# to retain these modifications.
task Funcotate {
    input {
        File ref_fasta
        File ref_fasta_index
        File ref_dict

        File input_vcf
        File input_vcf_idx

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
        Array[String]? annotation_defaults
        Array[String]? annotation_overrides
        Array[String]? funcotator_excluded_fields
        Boolean? filter_funcotations
        File? interval_list

        String? extra_args
        String? gcs_project_for_requester_pays

        # Runtime options:
        String gatk_docker

        File? gatk_override
        Int? mem
        Int? preemptible_attempts
        Int? max_retries
        Int? disk_space_gb
        Int? cpu

        Boolean use_ssd = true
    }

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int default_ram_mb = 1024 * 3
    Int machine_mem = if defined(mem) then mem * 1024 else default_ram_mb
    Int command_mem = machine_mem - 1024

    # Calculate disk size:
    Float ref_size_gb = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + size(ref_dict, "GiB")
    Float vcf_size_gb = size(input_vcf, "GiB") + size(input_vcf_idx, "GiB")
    Float datasources_size_gb = 25

    Int default_disk_space_gb = ceil(ref_size_gb + (datasources_size_gb * 2) + vcf_size_gb) + 20

    # Process input args:
    String output_maf = output_file_base_name + ".maf"
    String output_maf_index = output_maf + ".idx"
    String output_vcf = output_file_base_name + if compress then ".vcf.gz" else ".vcf"
    String output_vcf_idx = output_vcf +  if compress then ".tbi" else ".idx"
    String output_file = if output_format == "MAF" then output_maf else output_vcf
    String output_file_index = if output_format == "MAF" then output_maf_index else output_vcf_idx
    String transcript_selection_arg = if defined(transcript_selection_list) then " --transcript-list " else ""
    String annotation_def_arg = if defined(annotation_defaults) then " --annotation-default " else ""
    String annotation_over_arg = if defined(annotation_overrides) then " --annotation-override " else ""
    String filter_funcotations_args = if defined(filter_funcotations) && (filter_funcotations) then " --remove-filtered-variants " else ""
    String excluded_fields_args = if defined(funcotator_excluded_fields) then " --exclude-field " else ""
    String interval_list_arg = if defined(interval_list) then " -L " else ""
    String extra_args_arg = select_first([extra_args, ""])

    # Silly hack to allow us to use the dollar sign in the command section:
    String dollar = "$"

    parameter_meta{
        ref_fasta: { localization_optional: true }
        ref_fasta_index: { localization_optional: true }
        ref_dict: { localization_optional: true }
        input_vcf: { localization_optional: true }
        input_vcf_idx: { localization_optional: true }
    }

    command <<<
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

        # Hack to validate our WDL inputs:
        #
        # NOTE: This happens here so that we don't waste time copying down the data sources if there's an error.
        if [[ "~{output_format}" != "MAF" ]] && [[ "~{output_format}" != "VCF" ]] ; then
            echo "ERROR: Output format must be MAF or VCF."
        fi

        # download and extract Funcotator data sources
        gatk FuncotatorDataSourceDownloader \
            --somatic \
            --validate-integrity \
            --hg38 \
            --extract-after-download \
            --verbosity WARNING \
            ~{"--gcs-project-for-requester-pays " + gcs_project_for_requester_pays}

        # find where the datasources were extracted to
        DATA_SOURCES_FOLDER=$(find . -type d -name "funcotator_dataSources*" -print -quit)

        if ~{use_gnomad} ; then
            echo "Enabling gnomAD..."
            for potential_gnomad_gz in gnomAD_exome.tar.gz gnomAD_genome.tar.gz ; do
                if [[ -f ~{dollar}{DATA_SOURCES_FOLDER}/~{dollar}{potential_gnomad_gz} ]] ; then
                    cd ~{dollar}{DATA_SOURCES_FOLDER}
                    tar -zvxf ~{dollar}{potential_gnomad_gz}
                    cd -
                else
                    echo "ERROR: Cannot find gnomAD folder: ~{dollar}{potential_gnomad_gz}" 1>&2
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
            ~{annotation_def_arg}~{default="" sep=" --annotation-default " annotation_defaults} \
            ~{annotation_over_arg}~{default="" sep=" --annotation-override " annotation_overrides} \
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
        docker: gatk_docker
        bootDiskSizeGb: 12
        memory: machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space_gb, default_disk_space_gb]) + if use_ssd then " SSD" else " HDD"
        preemptible: select_first([preemptible_attempts, 3])
        maxRetries: select_first([max_retries, 0])
        cpu: select_first([cpu, 1])
    }

    output {
        File funcotated_output_file = "~{output_file}"
        File funcotated_output_file_index = "~{output_file_index}"
    }
}
