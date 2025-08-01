version 1.0

# Adapted from https://github.com/broadinstitute/gatk/blob/4.5.0.0/scripts/mutect2_wdl/mutect2.wdl
#
# Modifications:
#   - Remove tasks and options related to mutect3 training data
#   - Use SSDs instead of HDDs by default
#
# Use this file as a base and manually implement changes in the upstream code in order
# to retain these modifications.

workflow call_mutations {
    input {
        # basic inputs
        String sample_id
        File? intervals
        File ref_fasta
        File ref_fai
        File ref_dict
        File tumor_reads
        File tumor_reads_index
        File? normal_reads
        File? normal_reads_index

        # optional but usually recommended resources
        File? pon
        File? pon_idx
        File? gnomad
        File? gnomad_idx
        File? variants_for_contamination
        File? variants_for_contamination_idx

        # extra arguments
        String? m2_extra_args
        String? m2_extra_filtering_args
        String? getpileupsummaries_extra_args
        String? split_intervals_extra_args

        # additional modes and outputs
        File? realignment_index_bundle
        String? realignment_extra_args
        Boolean run_orientation_bias_mixture_model_filter = false
        Boolean make_bamout = false
        Boolean compress_vcfs = false
        File? gga_vcf
        File? gga_vcf_idx

        # runtime
        String gatk_docker
        File? gatk_override
        Int scatter_count
        Int preemptible = 1
        Int max_retries = 1
        Int small_task_cpu = 1
        Int small_task_mem = 4
        Int small_task_disk = 15
        Int learn_read_orientation_mem = 8000
        Int filter_alignment_artifacts_mem = 9000

        # Use as a last resort to increase the disk given to every task in case of ill behaving data
        Int emergency_extra_disk = 0
    }

    # Disk sizes used for dynamic sizing
    Int ref_size = ceil(size(ref_fasta, "GB") + size(ref_dict, "GB") + size(ref_fai, "GB"))
    Int tumor_reads_size = ceil(size(tumor_reads, "GB") + size(tumor_reads_index, "GB"))
    Int gnomad_vcf_size = if defined(gnomad) then ceil(size(gnomad, "GB")) else 0
    Int normal_reads_size = if defined(normal_reads) then ceil(size(normal_reads, "GB") + size(normal_reads_index, "GB")) else 0

    # This is added to every task as padding, should increase if systematically you need more disk for every call
    Int disk_pad = 10 + emergency_extra_disk

    Runtime standard_runtime = {"gatk_docker": gatk_docker, "gatk_override": gatk_override,
            "max_retries": max_retries, "preemptible": preemptible, "cpu": small_task_cpu,
            "machine_mem": small_task_mem * 1000, "command_mem": small_task_mem * 1000 - 500,
            "disk": small_task_disk + disk_pad}

    Int tumor_reads_size = ceil(size(tumor_reads, "GB") + size(tumor_reads_index, "GB"))
    Int normal_reads_size = if defined(normal_reads) then ceil(size(normal_reads, "GB") + size(normal_reads_index, "GB")) else 0

    Int m2_output_size = tumor_reads_size / scatter_count
    Int m2_per_scatter_size = ref_size + gnomad_vcf_size + m2_output_size + disk_pad

    call SplitIntervals {
        input:
            intervals = intervals,
            ref_fasta = ref_fasta,
            ref_fai = ref_fai,
            ref_dict = ref_dict,
            scatter_count = scatter_count,
            split_intervals_extra_args = split_intervals_extra_args,
            runtime_params = standard_runtime
    }

    scatter (subintervals in SplitIntervals.interval_files) {
        call M2 {
            input:
                intervals = subintervals,
                ref_fasta = ref_fasta,
                ref_fai = ref_fai,
                ref_dict = ref_dict,
                tumor_reads = tumor_reads,
                tumor_reads_index = tumor_reads_index,
                normal_reads = normal_reads,
                normal_reads_index = normal_reads_index,
                pon = pon,
                pon_idx = pon_idx,
                gnomad = gnomad,
                gnomad_idx = gnomad_idx,
                mem = small_task_mem,
                preemptible = 2,
                max_retries = max_retries,
                m2_extra_args = m2_extra_args,
                getpileupsummaries_extra_args = getpileupsummaries_extra_args,
                variants_for_contamination = variants_for_contamination,
                variants_for_contamination_idx = variants_for_contamination_idx,
                make_bamout = make_bamout,
                run_ob_filter = run_orientation_bias_mixture_model_filter,
                compress_vcfs = compress_vcfs,
                gga_vcf = gga_vcf,
                gga_vcf_idx = gga_vcf_idx,
                gatk_override = gatk_override,
                gatk_docker = gatk_docker,
                disk_space = m2_per_scatter_size
        }
    }

    Int merged_vcf_size = ceil(size(M2.unfiltered_vcf, "GB"))
    Int merged_bamout_size = ceil(size(M2.output_bamOut, "GB"))

    if (run_orientation_bias_mixture_model_filter) {
        call LearnReadOrientationModel {
            input:
                f1r2_tar_gz = M2.f1r2_counts,
                runtime_params = standard_runtime,
                mem = learn_read_orientation_mem
        }
    }

    call MergeVCFs {
        input:
            input_vcfs = M2.unfiltered_vcf,
            input_vcf_indices = M2.unfiltered_vcf_idx,
            compress_vcfs = compress_vcfs,
            runtime_params = standard_runtime
    }

    if (make_bamout) {
        call MergeBamOuts {
            input:
                ref_fasta = ref_fasta,
                ref_fai = ref_fai,
                ref_dict = ref_dict,
                bam_outs = M2.output_bamOut,
                runtime_params = standard_runtime,
                disk_space = ceil(merged_bamout_size * 4) + disk_pad,
        }
    }

    call MergeStats { input: stats = M2.stats, runtime_params = standard_runtime }

    if (defined(variants_for_contamination)) {
        call MergePileupSummaries as MergeTumorPileups {
            input:
                input_tables = flatten(M2.tumor_pileups),
                output_name = "tumor-pileups",
                ref_dict = ref_dict,
                runtime_params = standard_runtime
        }

        if (defined(normal_reads)){
            call MergePileupSummaries as MergeNormalPileups {
                input:
                    input_tables = flatten(M2.normal_pileups),
                    output_name = "normal-pileups",
                    ref_dict = ref_dict,
                    runtime_params = standard_runtime
            }
        }

        call CalculateContamination {
            input:
                tumor_pileups = MergeTumorPileups.merged_table,
                normal_pileups = MergeNormalPileups.merged_table,
                runtime_params = standard_runtime
        }
    }

    call Filter {
        input:
            sample_id = sample_id,
            ref_fasta = ref_fasta,
            ref_fai = ref_fai,
            ref_dict = ref_dict,
            intervals = intervals,
            unfiltered_vcf = MergeVCFs.merged_vcf,
            unfiltered_vcf_idx = MergeVCFs.merged_vcf_idx,
            compress_vcfs = compress_vcfs,
            mutect_stats = MergeStats.merged_stats,
            contamination_table = CalculateContamination.contamination_table,
            maf_segments = CalculateContamination.maf_segments,
            artifact_priors_tar_gz = LearnReadOrientationModel.artifact_prior_table,
            m2_extra_filtering_args = m2_extra_filtering_args,
            runtime_params = standard_runtime,
            disk_space = ceil(size(MergeVCFs.merged_vcf, "GB") * 4) + disk_pad
    }

    if (defined(realignment_index_bundle)) {
        call FilterAlignmentArtifacts {
            input:
                sample_id = sample_id,
                ref_fasta = ref_fasta,
                ref_fai = ref_fai,
                ref_dict = ref_dict,
                reads = tumor_reads,
                reads_index = tumor_reads_index,
                realignment_index_bundle = select_first([realignment_index_bundle]),
                realignment_extra_args = realignment_extra_args,
                compress_vcfs = compress_vcfs,
                input_vcf = Filter.filtered_vcf,
                input_vcf_idx = Filter.filtered_vcf_idx,
                runtime_params = standard_runtime,
                mem = filter_alignment_artifacts_mem
        }
    }

    output {
        File filtered_vcf = select_first([FilterAlignmentArtifacts.filtered_vcf, Filter.filtered_vcf])
        File filtered_vcf_idx = select_first([FilterAlignmentArtifacts.filtered_vcf_idx, Filter.filtered_vcf_idx])
        File filtering_stats = Filter.filtering_stats
    }
}

struct Runtime {
    String gatk_docker
    File? gatk_override
    Int max_retries
    Int preemptible
    Int cpu
    Int machine_mem
    Int command_mem
    Int disk
}

task SplitIntervals {
    input {
      File? intervals
      File ref_fasta
      File ref_fai
      File ref_dict
      Int scatter_count
      String? split_intervals_extra_args

      # runtime
      Runtime runtime_params
    }

    command <<<
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.gatk_override}

        mkdir interval-files
        java -Xmx~{runtime_params.command_mem}m -jar ${GATK_LOCAL_JAR} SplitIntervals \
            -R ~{ref_fasta} \
            ~{"-L " + intervals} \
            -scatter ~{scatter_count} \
            -O interval-files \
            ~{split_intervals_extra_args}
        cp interval-files/*.interval_list .
    >>>

    runtime {
        docker: runtime_params.gatk_docker
        memory: runtime_params.machine_mem + " MB"
        disks: "local-disk " + runtime_params.disk + " SSD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }

    output {
        Array[File] interval_files = glob("*.interval_list")
    }
}

task M2 {
    input {
        File? intervals
        File ref_fasta
        File ref_fai
        File ref_dict
        File tumor_reads
        File tumor_reads_index
        File? normal_reads
        File? normal_reads_index
        File? pon
        File? pon_idx
        File? gnomad
        File? gnomad_idx
        String? m2_extra_args
        String? getpileupsummaries_extra_args
        Boolean? make_bamout
        Boolean? run_ob_filter
        Boolean compress_vcfs
        File? gga_vcf
        File? gga_vcf_idx
        File? variants_for_contamination
        File? variants_for_contamination_idx

        File? gatk_override

        # runtime
        String gatk_docker
        Int? mem
        Int? preemptible
        Int? max_retries
        Int? disk_space
        Int? cpu
    }

    String output_vcf = "output" + if compress_vcfs then ".vcf.gz" else ".vcf"
    String output_vcf_idx = output_vcf + if compress_vcfs then ".tbi" else ".idx"

    String output_stats = output_vcf + ".stats"

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem) then mem * 1000 else 3500
    Int command_mem = machine_mem - 500

    parameter_meta{
        intervals: {localization_optional: true}
        ref_fasta: {localization_optional: true}
        ref_fai: {localization_optional: true}
        ref_dict: {localization_optional: true}
        tumor_reads: {localization_optional: true}
        tumor_reads_index: {localization_optional: true}
        normal_reads: {localization_optional: true}
        normal_reads_index: {localization_optional: true}
        pon: {localization_optional: true}
        pon_idx: {localization_optional: true}
        gnomad: {localization_optional: true}
        gnomad_idx: {localization_optional: true}
        gga_vcf: {localization_optional: true}
        gga_vcf_idx: {localization_optional: true}
        variants_for_contamination: {localization_optional: true}
        variants_for_contamination_idx: {localization_optional: true}
    }

    command <<<
        set -e

        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}
        export GCS_REQUESTER_PAYS_PROJECT="$(gcloud config get-value project -q)"

        # We need to create these files regardless, even if they stay empty
        touch bamout.bam
        touch f1r2.tar.gz
        touch dataset.txt
        echo "" > normal_name.txt

        java -Xmx~{command_mem}m -jar ${GATK_LOCAL_JAR} GetSampleName -R ~{ref_fasta} -I ~{tumor_reads} -O tumor_name.txt -encode \
        --gcs-project-for-requester-pays "${GCS_REQUESTER_PAYS_PROJECT}"
        tumor_command_line="-I ~{tumor_reads} -tumor `cat tumor_name.txt`"

        if [[ ! -z "~{normal_reads}" ]]; then
            java -Xmx~{command_mem}m -jar ${GATK_LOCAL_JAR} GetSampleName -R ~{ref_fasta} -I ~{normal_reads} -O normal_name.txt -encode \
            --gcs-project-for-requester-pays "${GCS_REQUESTER_PAYS_PROJECT}"
            normal_command_line="-I ~{normal_reads} -normal `cat normal_name.txt`"
        fi

        java -Xmx~{command_mem}m -jar ${GATK_LOCAL_JAR} Mutect2 \
            -R ~{ref_fasta} \
            $tumor_command_line \
            $normal_command_line \
            ~{"--germline-resource " + gnomad} \
            ~{"-pon " + pon} \
            ~{"-L " + intervals} \
            ~{"--alleles " + gga_vcf} \
            -O "~{output_vcf}" \
            ~{true='--bam-output bamout.bam' false='' make_bamout} \
            ~{true='--f1r2-tar-gz f1r2.tar.gz' false='' run_ob_filter} \
            ~{m2_extra_args} \
            --gcs-project-for-requester-pays "${GCS_REQUESTER_PAYS_PROJECT}"

        m2_exit_code=$?

        ### GetPileupSummaries

        # If the variants for contamination and the intervals for this scatter don't intersect, GetPileupSummaries
        # throws an error.  However, there is nothing wrong with an empty intersection for our purposes; it simply doesn't
        # contribute to the merged pileup summaries that we create downstream.  We implement this by with array outputs.
        # If the tool errors, no table is created and the glob yields an empty array.
        set +e

        if [[ ! -z "~{variants_for_contamination}" ]]; then
            java -Xmx~{command_mem}m -jar ${GATK_LOCAL_JAR} GetPileupSummaries -R ~{ref_fasta} -I ~{tumor_reads} ~{"--interval-set-rule INTERSECTION -L " + intervals} \
                -V ~{variants_for_contamination} -L ~{variants_for_contamination} -O tumor-pileups.table ~{getpileupsummaries_extra_args} \
                --gcs-project-for-requester-pays "${GCS_REQUESTER_PAYS_PROJECT}"


            if [[ ! -z "~{normal_reads}" ]]; then
                java -Xmx~{command_mem}m -jar ${GATK_LOCAL_JAR} GetPileupSummaries -R ~{ref_fasta} -I ~{normal_reads} ~{"--interval-set-rule INTERSECTION -L " + intervals} \
                    -V ~{variants_for_contamination} -L ~{variants_for_contamination} -O normal-pileups.table ~{getpileupsummaries_extra_args} \
                    --gcs-project-for-requester-pays "${GCS_REQUESTER_PAYS_PROJECT}"
            fi
        fi

        # the script only fails if Mutect2 itself fails
        exit $m2_exit_code
    >>>

    runtime {
        docker: gatk_docker
        memory: machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space, 100]) + " SSD"
        preemptible: select_first([preemptible, 10])
        maxRetries: select_first([max_retries, 0])
        cpu: select_first([cpu, 1])
    }

    output {
        File unfiltered_vcf = "~{output_vcf}"
        File unfiltered_vcf_idx = "~{output_vcf_idx}"
        File output_bamOut = "bamout.bam"
        String tumor_sample = read_string("tumor_name.txt")
        String normal_sample = read_string("normal_name.txt")
        File stats = "~{output_stats}"
        File f1r2_counts = "f1r2.tar.gz"
        Array[File] tumor_pileups = glob("*tumor-pileups.table")
        Array[File] normal_pileups = glob("*normal-pileups.table")
    }
}

task MergeVCFs {
    input {
      Array[File] input_vcfs
      Array[File] input_vcf_indices
      Boolean compress_vcfs
      Runtime runtime_params
    }

    String output_vcf = if compress_vcfs then "merged.vcf.gz" else "merged.vcf"
    String output_vcf_idx = output_vcf + if compress_vcfs then ".tbi" else ".idx"

    # using MergeVcfs instead of GatherVcfs so we can create indices
    # WARNING 2015-10-28 15:01:48 GatherVcfs  Index creation not currently supported when gathering block compressed VCFs.
    command <<<
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.gatk_override}
        java -Xmx~{runtime_params.command_mem}m -jar ${GATK_LOCAL_JAR} MergeVcfs -I ~{sep=' -I ' input_vcfs} -O ~{output_vcf}
    >>>

    runtime {
        docker: runtime_params.gatk_docker
        memory: runtime_params.machine_mem + " MB"
        disks: "local-disk " + runtime_params.disk + " SSD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }

    output {
        File merged_vcf = "~{output_vcf}"
        File merged_vcf_idx = "~{output_vcf_idx}"
    }
}

task MergeBamOuts {
    input {
      File ref_fasta
      File ref_fai
      File ref_dict
      Array[File]+ bam_outs
      Runtime runtime_params
      Int? disk_space   #override to request more disk than default small task params
    }

    command <<<
        # This command block assumes that there is at least one file in bam_outs.
        #  Do not call this task if len(bam_outs) == 0
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.gatk_override}
        java -Xmx~{runtime_params.command_mem}m -jar ${GATK_LOCAL_JAR} GatherBamFiles \
            -I ~{sep=" -I " bam_outs} -O unsorted.out.bam -R ~{ref_fasta}

        # We must sort because adjacent scatters may have overlapping (padded) assembly regions, hence
        # overlapping bamouts

        java -Xmx~{runtime_params.command_mem}m -jar ${GATK_LOCAL_JAR} SortSam -I unsorted.out.bam \
            -O bamout.bam --SORT_ORDER coordinate -VALIDATION_STRINGENCY LENIENT
        java -Xmx~{runtime_params.command_mem}m -jar ${GATK_LOCAL_JAR} BuildBamIndex -I bamout.bam -VALIDATION_STRINGENCY LENIENT
    >>>

    runtime {
        docker: runtime_params.gatk_docker
        memory: runtime_params.machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space, runtime_params.disk]) + " SSD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }

    output {
        File merged_bam_out = "bamout.bam"
        File merged_bam_out_index = "bamout.bai"
    }
}


task MergeStats {
    input {
      Array[File]+ stats
      Runtime runtime_params
    }

    command <<<
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.gatk_override}

        java -Xmx~{runtime_params.command_mem}m -jar ${GATK_LOCAL_JAR} MergeMutectStats \
            -stats ~{sep=" -stats " stats} -O merged.stats
    >>>

    runtime {
        docker: runtime_params.gatk_docker
        memory: runtime_params.machine_mem + " MB"
        disks: "local-disk " + runtime_params.disk + " SSD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }

    output {
        File merged_stats = "merged.stats"
    }
}

task MergePileupSummaries {
    input {
      Array[File] input_tables
      String output_name
      File ref_dict
      Runtime runtime_params
    }

    command <<<
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.gatk_override}

        java -Xmx~{runtime_params.command_mem}m -jar ${GATK_LOCAL_JAR} GatherPileupSummaries \
        --sequence-dictionary ~{ref_dict} \
        -I ~{sep=' -I ' input_tables} \
        -O ~{output_name}.tsv
    >>>

    runtime {
        docker: runtime_params.gatk_docker
        memory: runtime_params.machine_mem + " MB"
        disks: "local-disk " + runtime_params.disk + " SSD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }

    output {
        File merged_table = "~{output_name}.tsv"
    }
}

# Learning step of the orientation bias mixture model, which is the recommended orientation bias filter as of September 2018
task LearnReadOrientationModel {
    input {
      Array[File] f1r2_tar_gz
      Runtime runtime_params
      Int? mem  #override memory
    }

    Int machine_mem = select_first([mem, runtime_params.machine_mem])
    Int command_mem = machine_mem - 1000

    command <<<
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.gatk_override}

        java -Xmx~{command_mem}m -jar ${GATK_LOCAL_JAR} LearnReadOrientationModel \
            -I ~{sep=" -I " f1r2_tar_gz} \
            -O "artifact-priors.tar.gz"
    >>>

    runtime {
        docker: runtime_params.gatk_docker
        memory: machine_mem + " MB"
        disks: "local-disk " + runtime_params.disk + " SSD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }

    output {
        File artifact_prior_table = "artifact-priors.tar.gz"
    }

}

task CalculateContamination {
    input {
      String? intervals
      File tumor_pileups
      File? normal_pileups
      Runtime runtime_params
    }

    command <<<
        set -e

        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.gatk_override}

        java -Xmx~{runtime_params.command_mem}m -jar ${GATK_LOCAL_JAR} CalculateContamination -I ~{tumor_pileups} \
        -O contamination.table --tumor-segmentation segments.table ~{"-matched " + normal_pileups}
    >>>

    runtime {
        docker: runtime_params.gatk_docker
        memory: runtime_params.machine_mem + " MB"
        disks: "local-disk " + runtime_params.disk + " SSD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }

    output {
        File contamination_table = "contamination.table"
        File maf_segments = "segments.table"
    }
}

task Filter {
    input {
        String sample_id
        File? intervals
        File ref_fasta
        File ref_fai
        File ref_dict
        File unfiltered_vcf
        File unfiltered_vcf_idx
        Boolean compress_vcfs
        File? mutect_stats
        File? artifact_priors_tar_gz
        File? contamination_table
        File? maf_segments
        String? m2_extra_filtering_args

        Runtime runtime_params
        Int? disk_space
    }

    String output_vcf = if compress_vcfs then "~{sample_id}.filtered.vcf.gz" else "~{sample_id}.filtered.vcf"
    String output_vcf_idx = output_vcf + if compress_vcfs then ".tbi" else ".idx"

    parameter_meta{
      ref_fasta: {localization_optional: true}
      ref_fai: {localization_optional: true}
      ref_dict: {localization_optional: true}
    }

    command <<<
        set -e

        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.gatk_override}

        java -Xmx~{runtime_params.command_mem}m -jar ${GATK_LOCAL_JAR} FilterMutectCalls -V ~{unfiltered_vcf} \
            -R ~{ref_fasta} \
            -O "~{output_vcf}" \
            ~{"--contamination-table " + contamination_table} \
            ~{"--tumor-segmentation " + maf_segments} \
            ~{"--ob-priors " + artifact_priors_tar_gz} \
            ~{"-stats " + mutect_stats} \
            --filtering-stats "~{sample_id}.filtering.stats" \
            ~{m2_extra_filtering_args}
    >>>

    runtime {
        docker: runtime_params.gatk_docker
        memory: runtime_params.machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space, runtime_params.disk]) + " SSD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }

    output {
        File filtered_vcf = "~{output_vcf}"
        File filtered_vcf_idx = "~{output_vcf_idx}"
        File filtering_stats = "~{sample_id}.filtering.stats"
    }
}

task FilterAlignmentArtifacts {
    input {
        String sample_id
        File ref_fasta
        File ref_fai
        File ref_dict
        File input_vcf
        File input_vcf_idx
        File reads
        File reads_index
        Boolean compress_vcfs
        File realignment_index_bundle
        String? realignment_extra_args
        Runtime runtime_params
        Int mem
    }

    String output_vcf = if compress_vcfs then "~{sample_id}.filtered.vcf.gz" else "~{sample_id}.filtered.vcf"
    String output_vcf_idx = output_vcf +  if compress_vcfs then ".tbi" else ".idx"

    Int machine_mem = mem
    Int command_mem = machine_mem - 500

    parameter_meta{
      ref_fasta: {localization_optional: true}
      ref_fai: {localization_optional: true}
      ref_dict: {localization_optional: true}
      input_vcf: {localization_optional: true}
      input_vcf_idx: {localization_optional: true}
      reads: {localization_optional: true}
      reads_index: {localization_optional: true}
    }

    command <<<
        set -e

        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.gatk_override}
        export GCS_REQUESTER_PAYS_PROJECT="$(gcloud config get-value project -q)"

        java -Xmx~{command_mem}m -jar ${GATK_LOCAL_JAR} FilterAlignmentArtifacts \
            -R ~{ref_fasta} \
            -V ~{input_vcf} \
            -I ~{reads} \
            --bwa-mem-index-image ~{realignment_index_bundle} \
            ~{realignment_extra_args} \
            -O ~{output_vcf} \
            --gcs-project-for-requester-pays "${GCS_REQUESTER_PAYS_PROJECT}"
    >>>

    runtime {
        docker: runtime_params.gatk_docker
        memory: machine_mem + " MB"
        disks: "local-disk " + runtime_params.disk + " SSD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }

    output {
        File filtered_vcf = "~{output_vcf}"
        File filtered_vcf_idx = "~{output_vcf_idx}"
    }
}
