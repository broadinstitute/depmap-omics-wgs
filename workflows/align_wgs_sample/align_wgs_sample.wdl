version 1.0

# Adapted from https://github.com/broadinstitute/warp/blob/GDCWholeGenomeSomaticSingleSample_v1.3.1/pipelines/broad/dna_seq/somatic/single_sample/wgs/gdc_genome/GDCWholeGenomeSomaticSingleSample.wdl
#
# Modifications:
#   - Separate CRAM/CRAI vs. BAM/BAI workflow inputs consolidated into singular
#     `input_cram_bam` and `input_crai_bai` ones, with necessary logic to handle either
#     file type.
#   - Use SSDs instead of HDDs for faster (de)localization and higher compute instance
#     availability.
#   - No `ubam` input option.
#   - `CheckContamination` task copied into this file.
#   - `CheckContamination` and `collect_insert_size_metrics` made optional.
#   - Remove `validation_report` and `unmapped_bams` workflow outputs.
#   - Scatter BQSR tasks.
#   - Check to see if BQSR has already been run but the input BAM file is missing `OQ`
#     tags, which would make re-running BQSR inappropriate.
#   - Added optional ref_alt input for BWA.
#
# Use this file as a base and manually implement changes in the upstream code in order
# to retain these modifications.

import "./include/CramToUnmappedBams.wdl" as ToUbams

struct FastqPairRecord {
    File forward_fastq
    File reverse_fastq
    String readgroup
    String readgroup_id
}

struct FastqSingleRecord {
    File fastq
    String readgroup
    String readgroup_id
}

workflow align_wgs_sample {
    input {
        String sample_id
        String input_type # "CRAM", "BAM", or "FASTQ"
        File? input_cram_bam
        File? input_crai_bai
        Array[File]? fastqs_pe_forward
        Array[File]? fastqs_pe_reverse
        Array[String]? fastqs_pe_readgroup
        Array[String]? fastqs_pe_readgroup_id
        Array[File]? fastqs_se
        Array[String]? fastqs_se_readgroup
        Array[String]? fastqs_se_readgroup_id
        File? delivery_ref_fasta
        File? delivery_ref_fasta_index
        File? output_map
        String? unmapped_bam_suffix

        Boolean do_check_contamination = false
        Boolean do_collect_insert_size_metrics = false
        File? contamination_vcf
        File? contamination_vcf_index
        Boolean run_bqsr = true
        Boolean rerun_bqsr = true
        File dbsnp_vcf
        File dbsnp_vcf_index

        File ref_fasta
        File ref_fasta_index
        File ref_dict
        File ref_amb
        File ref_ann
        File ref_bwt
        File ref_pac
        File ref_sa
        File? ref_alt
    }

    if (input_type == "BAM" || input_type == "CRAM") {
        call ToUbams.CramToUnmappedBams {
            input:
                input_type = input_type,
                input_cram_bam = select_first([input_cram_bam]),
                input_crai_bai = select_first([input_crai_bai]),
                ref_fasta = select_first([delivery_ref_fasta, ref_fasta]),
                ref_fasta_index = select_first([delivery_ref_fasta_index, ref_fasta_index]),
                output_map = output_map,
                unmapped_bam_suffix = unmapped_bam_suffix
        }

        scatter (ubam in CramToUnmappedBams.unmapped_bams) {
            call bam_readgroup_to_contents {
                input: bam = ubam
            }

            call biobambam_bamtofastq {
                input: filename = ubam
            }
        }

        Array[Object] readgroups = flatten(bam_readgroup_to_contents.readgroups)

        Array[File] fastq1 = flatten(biobambam_bamtofastq.output_fastq1)
        Array[File] fastq2 = flatten(biobambam_bamtofastq.output_fastq2)
        Array[File] fastq_o1 = flatten(biobambam_bamtofastq.output_fastq_o1)
        Array[File] fastq_o2 = flatten(biobambam_bamtofastq.output_fastq_o2)
        Array[File] fastq_s = flatten(biobambam_bamtofastq.output_fastq_s)

        Int pe_count = length(fastq1)
        Int o1_count = length(fastq_o1)
        Int o2_count = length(fastq_o2)
        Int s_count = length(fastq_s)

        if (pe_count > 0) {
            call emit_pe_records {
                input:
                    fastq1_files = fastq1,
                    fastq2_files = fastq2,
                    readgroups = readgroups
            }

            scatter (pe_record in emit_pe_records.fastq_pair_records) {
                call bwa_pe {
                    input:
                        fastq_record = pe_record,
                        ref_fasta = ref_fasta,
                        ref_fasta_index = ref_fasta_index,
                        ref_dict = ref_dict,
                        ref_amb = ref_amb,
                        ref_ann = ref_ann,
                        ref_bwt = ref_bwt,
                        ref_pac = ref_pac,
                        ref_sa = ref_sa,
                        ref_alt = ref_alt
                }
            }
        }

        if (o1_count + o2_count + s_count > 0) {
            call emit_se_records {
                input:
                    fastq_o1_files = fastq_o1,
                    fastq_o2_files = fastq_o2,
                    fastq_s_files = fastq_s,
                    readgroups = readgroups
            }

            scatter (se_record in emit_se_records.fastq_single_records) {
                call bwa_se {
                    input:
                        fastq_record = se_record,
                        ref_fasta = ref_fasta,
                        ref_fasta_index = ref_fasta_index,
                        ref_dict = ref_dict,
                        ref_amb = ref_amb,
                        ref_ann = ref_ann,
                        ref_bwt = ref_bwt,
                        ref_pac = ref_pac,
                        ref_sa = ref_sa,
                        ref_alt = ref_alt
                }
            }
        }
    }

    if (input_type == "FASTQ") {
        Int fastqs_pe_count = length(select_first([fastqs_pe_forward]))

        if (fastqs_pe_count > 0) {
            call emit_fastq_pe_records {
                input:
                    fastqs_pe_forward = select_first([fastqs_pe_forward]),
                    fastqs_pe_reverse = select_first([fastqs_pe_reverse]),
                    fastqs_pe_readgroup = select_first([fastqs_pe_readgroup]),
                    fastqs_pe_readgroup_id = select_first([fastqs_pe_readgroup_id])
            }

            scatter (pe_record in emit_fastq_pe_records.fastq_pair_records) {
                call bwa_pe as fastq_bwa_pe {
                    input:
                        fastq_record = pe_record,
                        ref_fasta = ref_fasta,
                        ref_fasta_index = ref_fasta_index,
                        ref_dict = ref_dict,
                        ref_amb = ref_amb,
                        ref_ann = ref_ann,
                        ref_bwt = ref_bwt,
                        ref_pac = ref_pac,
                        ref_sa = ref_sa,
                        ref_alt = ref_alt
                }
            }
        }

        Int fastqs_se_count = length(select_first([fastqs_se]))

        if (fastqs_se_count > 0) {
            call emit_fastqs_se_records {
                input:
                    fastqs_se = select_first([fastqs_se]),
                    fastqs_se_readgroup = select_first([fastqs_se_readgroup]),
                    fastqs_se_readgroup_id = select_first([fastqs_se_readgroup_id])
            }

            scatter (se_record in emit_fastqs_se_records.fastq_single_records) {
                call bwa_se as fastq_bwa_se {
                    input:
                        fastq_record = se_record,
                        ref_fasta = ref_fasta,
                        ref_fasta_index = ref_fasta_index,
                        ref_dict = ref_dict,
                        ref_amb = ref_amb,
                        ref_ann = ref_ann,
                        ref_bwt = ref_bwt,
                        ref_pac = ref_pac,
                        ref_sa = ref_sa,
                        ref_alt = ref_alt
                }
            }
        }
    }

    Array[File] aligned_bams = flatten([
        select_first([bwa_pe.bam, []]),
        select_first([bwa_se.bam, []]),
        select_first([fastq_bwa_pe.bam, []]),
        select_first([fastq_bwa_se.bam, []]),
    ])

    call picard_markduplicates {
        input:
            bams = aligned_bams,
            outbam = sample_id + ".mrkdup.bam"
    }

    call sort_and_index_markdup_bam {
        input:
            input_bam = picard_markduplicates.bam,
            output_bam_basename = sample_id + ".analysis_ready"
    }

    if (do_check_contamination) {
        call CalculateSomaticContamination as check_contamination {
            input:
                reference = ref_fasta,
                reference_dict = ref_dict,
                reference_index = ref_fasta_index,
                tumor_cram_or_bam = sort_and_index_markdup_bam.output_bam,
                tumor_crai_or_bai = sort_and_index_markdup_bam.output_bai,
                contamination_vcf = contamination_vcf,
                contamination_vcf_index = contamination_vcf_index
        }
    }

    if (run_bqsr) {
        if (input_type == "BAM" || input_type == "CRAM") {
            call check_bqsr_and_oq_tags {
                input:
                    bam = select_first([input_cram_bam]),
                    bai = select_first([input_crai_bai]),
                    ref_fasta = ref_fasta
            }
        }

        if (input_type == "FASTQ" || !check_bqsr_and_oq_tags.bqsr_performed || (rerun_bqsr && check_bqsr_and_oq_tags.has_oq_tags)) {
            call CreateSequenceGroupingTSV {
                input:
                    ref_dict = ref_dict
            }

            Int n_bqsr_splits = length(CreateSequenceGroupingTSV.sequence_grouping_with_unmapped)

            scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping) {
                call BaseRecalibrator {
                    input:
                        bqsr_regions = subgroup,
                        n_bqsr_splits = n_bqsr_splits,
                        bam = sort_and_index_markdup_bam.output_bam,
                        bam_index = sort_and_index_markdup_bam.output_bai,
                        ref_fasta = ref_fasta,
                        ref_fasta_index = ref_fasta_index,
                        ref_dict = ref_dict,
                        dbsnp_vcf = dbsnp_vcf,
                        dbsnp_vcf_index = dbsnp_vcf_index
                }
            }

            call GatherBqsrReports {
                input:
                    input_bqsr_reports = BaseRecalibrator.bqsr_recal_file,
                    output_report_filename = sample_id + ".bqsr.grp"
            }

            scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping_with_unmapped) {
                call ApplyBQSR {
                    input:
                        bqsr_regions = subgroup,
                        n_bqsr_splits = n_bqsr_splits,
                        bam = sort_and_index_markdup_bam.output_bam,
                        bam_index = sort_and_index_markdup_bam.output_bai,
                        bqsr_recal_file = GatherBqsrReports.output_bqsr_report
                }
            }

            call GatherBamFiles as BqsrGatherBamFiles {
                input:
                    input_bams = ApplyBQSR.recalibrated_bam,
                    output_bam_basename = sample_id + ".analysis_ready"
            }
        }
    }

    File final_bam = select_first([BqsrGatherBamFiles.output_bam, sort_and_index_markdup_bam.output_bam])
    File final_bai = select_first([BqsrGatherBamFiles.output_bai, sort_and_index_markdup_bam.output_bai])

    if (do_collect_insert_size_metrics) {
        call collect_insert_size_metrics {
            input:
                input_bam = final_bam,
                output_bam_basename = sample_id
        }
    }

    output {
        File analysis_ready_bam = final_bam
        File analysis_ready_bai = final_bai
        File md_metrics = picard_markduplicates.metrics
        File? insert_size_metrics = collect_insert_size_metrics.insert_size_metrics
        File? insert_size_histogram_pdf = collect_insert_size_metrics.insert_size_histogram_pdf
        File? contamination = check_contamination.contamination
    }

    meta {
        allowNestedInputs: true
    }
}

task bam_readgroup_to_contents {
    input {
        File bam
        Int preemptible = 5
        Int max_retries = 1
        Int cpu = 1
        Int additional_memory_mb = 0
        Int additional_disk_gb = 0
    }

    Int mem = ceil(size(bam, "MiB")) + 2000 + additional_memory_mb
    Int disk_space = ceil(size(bam, "GiB")) + 10 + additional_disk_gb

    parameter_meta {
        bam: { localization_optional: true }
    }

    command <<<
        set -euo pipefail

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

        samtools view -H ~{bam} | \
            awk 'BEGIN {                                                                \
                    OFS = FS = "\t";                                                    \
                    header = "ID\tBC\tCN\tDS\tDT\tFO\tKS\tLB\tPG\tPI\tPL\tPM\tPU\tSM";  \
                    split(header, header_ary, "\t");                                    \
                    for (i=1; i in header_ary; i++) {                                   \
                        header_pos[header_ary[i]] = i                                   \
                    };                                                                  \
                    print header                                                        \
                 }                                                                      \
                 /^@RG/ {                                                               \
                    for (i=2; i<=NF; i++) {                                             \
                        split($i, rg, ":");                                             \
                        row_ary[header_pos[rg[1]]] = rg[2];                             \
                    };                                                                  \
                    row = row_ary[1];                                                   \
                    for (i=2; i in header_ary; i++) {                                   \
                        row = row "\t";                                                 \
                        if (i in row_ary)                                               \
                            row = row row_ary[i];                                       \
                    };                                                                  \
                    delete row_ary;                                                     \
                    print row                                                           \
                 }'
    >>>

    output {
        Array[Object] readgroups = read_objects(stdout())
    }

    runtime {
        docker: "broadgdac/samtools:1.10"
        memory: mem + " MiB"
        disks: "local-disk " + disk_space + " SSD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }
}

task biobambam_bamtofastq {
    input {
        Int collate = 1
        String exclude = "QCFAIL,SECONDARY,SUPPLEMENTARY"
        File filename
        Int gz = 1
        String inputformat = "bam"
        Int level = 5
        String outputdir = "."
        Int outputperreadgroup = 1
        String outputperreadgroupsuffixF = "_1.fq.gz"
        String outputperreadgroupsuffixF2 = "_2.fq.gz"
        String outputperreadgroupsuffixO = "_o1.fq.gz"
        String outputperreadgroupsuffixO2 = "_o2.fq.gz"
        String outputperreadgroupsuffixS = "_s.fq.gz"
        Int tryoq = 1
        String T = "tempfq"
        Int cpu = 1
        Int preemptible = 2
        Int max_retries = 1
        Int additional_memory_mb = 0
        Int additional_disk_gb = 0
    }

    Int mem = ceil(size(filename, "MiB")) + 2000 + additional_memory_mb
    Int disk_space = ceil(size(filename, "GiB") * 2) + 10 + additional_disk_gb

    command <<<
        set -euo pipefail

        /usr/local/bin/bamtofastq \
            T=~{T} \
            collate=~{collate} \
            exclude=~{exclude} \
            filename=~{filename} \
            gz=~{gz} \
            inputformat=~{inputformat} \
            level=~{level} \
            outputdir=~{outputdir} \
            outputperreadgroup=~{outputperreadgroup} \
            outputperreadgroupsuffixF=~{outputperreadgroupsuffixF} \
            outputperreadgroupsuffixF2=~{outputperreadgroupsuffixF2} \
            outputperreadgroupsuffixO=~{outputperreadgroupsuffixO} \
            outputperreadgroupsuffixO2=~{outputperreadgroupsuffixO2} \
            outputperreadgroupsuffixS=~{outputperreadgroupsuffixS} \
            tryoq=~{tryoq}
    >>>

    output {
        Array[File] output_fastq1 = glob("*~{outputperreadgroupsuffixF}")
        Array[File] output_fastq2 = glob("*~{outputperreadgroupsuffixF2}")
        Array[File] output_fastq_o1 = glob("*~{outputperreadgroupsuffixO}")
        Array[File] output_fastq_o2 = glob("*~{outputperreadgroupsuffixO2}")
        Array[File] output_fastq_s = glob("*~{outputperreadgroupsuffixS}")
    }

    runtime {
        docker: "broadgdac/biobambam2:2.0.87-release-20180301132713"
        memory: mem + " MiB"
        disks: "local-disk " + disk_space + " SSD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }
}

task emit_pe_records {
    input {
        Array[File]+ fastq1_files
        Array[File]+ fastq2_files
        Array[Object]+ readgroups
        String fastq1_suffix = "_1.fq.gz"
        Int preemptible = 5
        Int max_retries = 1
        Int cpu = 1
        Int additional_memory_mb = 0
        Int additional_disk_gb = 0
    }

    Int mem = ceil(size(fastq1_files, "MiB") + size(fastq2_files, "MiB")) + 2000 + additional_memory_mb
    Int disk_space = ceil(size(fastq1_files, "GiB") + size(fastq2_files, "GiB")) + 10 + additional_disk_gb

    File readgroups_tsv = write_objects(readgroups)

    parameter_meta {
        fastq1_files: { localization_optional: true }
        fastq2_files: { localization_optional: true }
    }

    command <<<
        set -euo pipefail

        python <<CODE
        from csv import DictReader

        def basename(bucket_path):
            return bucket_path.rsplit('/', 1)[-1]

        fastq1_files = sorted("~{sep=',' fastq1_files}".split(','), key=basename)
        fastq2_files = sorted("~{sep=',' fastq2_files}".split(','), key=basename)

        readgroups = dict()
        with open("~{readgroups_tsv}") as readgroups_tsv:
            for readgroup in DictReader(readgroups_tsv, dialect="excel-tab"):
                readgroups[readgroup["ID"]] = r"@RG\t" + r"\t".join(
                    "{}:{}".format(key, value)
                    for key, value in readgroup.items() if value)

        print("forward_fastq\treverse_fastq\treadgroup\treadgroup_id")
        for fastq1, fastq2 in zip(fastq1_files, fastq2_files):
            rg_id = basename(fastq1).rsplit("~{fastq1_suffix}", 1)[0]
            print("\t".join([fastq1, fastq2, readgroups[rg_id], rg_id]))
        CODE
    >>>

    output {
        Array[FastqPairRecord] fastq_pair_records = read_objects(stdout())
    }

    runtime {
        docker: "python:3.8-slim"
        memory: mem + " MiB"
        disks: "local-disk " + disk_space + " SSD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }
}

task emit_se_records {
    input {
        Array[File] fastq_o1_files
        Array[File] fastq_o2_files
        Array[File] fastq_s_files
        Array[Object]+ readgroups
        String fastq_o1_suffix = "_o1.fq.gz"
        String fastq_o2_suffix = "_o2.fq.gz"
        String fastq_s_suffix = "_s.fq.gz"
        Int preemptible = 5
        Int max_retries = 1
        Int cpu = 1
        Int additional_memory_mb = 0
        Int additional_disk_gb = 0
    }

    Int mem = ceil(size(fastq_o1_files, "MiB") + size(fastq_o2_files, "MiB") + size(fastq_s_files, "MiB")) + 2000 + additional_memory_mb
    Int disk_space = ceil(size(fastq_o1_files, "GiB") + size(fastq_o2_files, "GiB") + size(fastq_s_files, "GiB")) + 10 + additional_disk_gb

    File readgroups_tsv = write_objects(readgroups)

    parameter_meta {
        fastq_o1_files: { localization_optional: true }
        fastq_o2_files: { localization_optional: true }
        fastq_s_files: { localization_optional: true }
    }

    command <<<
        set -euo pipefail

        python <<CODE
        from csv import DictReader

        def basename(bucket_path):
            return bucket_path.rsplit('/', 1)[-1]

        def emit_records(fastqs, suffix, readgroups):
            for fastq in fastqs:
                if not fastq:
                    return
                rg_id = basename(fastq).rsplit(suffix, 1)[0]
                print("\t".join([fastq, readgroups[rg_id], rg_id]))

        fastq_o1_files = "~{sep=',' fastq_o1_files}".split(',')
        fastq_o2_files = "~{sep=',' fastq_o2_files}".split(',')
        fastq_s_files = "~{sep=',' fastq_s_files}".split(',')

        readgroups = dict()
        with open("~{readgroups_tsv}") as readgroups_tsv:
            for readgroup in DictReader(readgroups_tsv, dialect="excel-tab"):
                readgroups[readgroup["ID"]] = r"@RG\t" + r"\t".join(
                    "{}:{}".format(key, value)
                    for key, value in readgroup.items() if value)

        print("fastq\treadgroup\treadgroup_id")
        emit_records(fastq_o1_files, "~{fastq_o1_suffix}", readgroups)
        emit_records(fastq_o2_files, "~{fastq_o2_suffix}", readgroups)
        emit_records(fastq_s_files, "~{fastq_s_suffix}", readgroups)
        CODE
    >>>

    output {
        Array[FastqSingleRecord] fastq_single_records = read_objects(stdout())
    }

    runtime {
        docker: "python:3.8-slim"
        memory: mem + " MiB"
        disks: "local-disk " + disk_space + " SSD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }
}

task emit_fastq_pe_records {
    input {
        Array[File] fastqs_pe_forward
        Array[File] fastqs_pe_reverse
        Array[String] fastqs_pe_readgroup
        Array[String] fastqs_pe_readgroup_id

        # any linux image works here
        String docker_image = "us-central1-docker.pkg.dev/depmap-omics/terra-images/bcftools"
        String docker_image_hash_or_tag = ":production"
        Int mem_gb = 2
        Int cpu = 1
        Int preemptible = 2
        Int max_retries = 1
        Int additional_disk_gb = 0
        Int disk_space = 10 # GB
    }

    parameter_meta {
        fastqs_pe_forward: { localization_optional: true }
        fastqs_pe_reverse: { localization_optional: true }
    }

    command <<<
        set -euo pipefail

        # Convert the WDL arrays (space-separated strings) into bash arrays.
        IFS=$' ' read -r -a fwd_arr  <<< "~{sep=' ' fastqs_pe_forward}"
        IFS=$' ' read -r -a rev_arr  <<< "~{sep=' ' fastqs_pe_reverse}"
        IFS=$' ' read -r -a rg_arr   <<< "~{sep=' ' fastqs_pe_readgroup}"
        IFS=$' ' read -r -a rgid_arr <<< "~{sep=' ' fastqs_pe_readgroup_id}"

        # Validate lengths
        len=${#fwd_arr[@]}
        if [[ ${#rev_arr[@]} -ne $len || ${#rg_arr[@]} -ne $len || ${#rgid_arr[@]} -ne $len ]]; then
            echo "ERROR: input array lengths differ" 1>&2
            exit 1
        fi

        # Write header
        printf "forward_fastq\treverse_fastq\treadgroup\treadgroup_id\n" > records.tsv

        # Emit rows
        for ((i=0; i<len; i++)); do
            printf "%s\t%s\t%s\t%s\n" \
                "${fwd_arr[i]}" \
                "${rev_arr[i]}" \
                "${rg_arr[i]}" \
                "${rgid_arr[i]}" >> records.tsv
        done
    >>>

    output {
        Array[FastqPairRecord] fastq_pair_records = read_objects("records.tsv")
    }

    runtime {
        docker: "~{docker_image}~{docker_image_hash_or_tag}"
        memory: "~{mem_gb} GB"
        disks: "local-disk ~{disk_space} SSD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }
}

task emit_fastqs_se_records {
    input {
        Array[File] fastqs_se
        Array[String] fastqs_se_readgroup
        Array[String] fastqs_se_readgroup_id

        # any linux image works here
        String docker_image = "us-central1-docker.pkg.dev/depmap-omics/terra-images/bcftools"
        String docker_image_hash_or_tag = ":production"
        Int mem_gb = 2
        Int cpu = 1
        Int preemptible = 2
        Int max_retries = 1
        Int additional_disk_gb = 0
        Int disk_space = 10 # GB
    }

    parameter_meta {
        fastqs_se: { localization_optional: true }
    }

    command <<<
        set -euo pipefail

        # Convert the WDL arrays (space-separated strings) into bash arrays.
        IFS=$' ' read -r -a fq_arr  <<< "~{sep=' ' fastqs_se}"
        IFS=$' ' read -r -a rg_arr   <<< "~{sep=' ' fastqs_se_readgroup}"
        IFS=$' ' read -r -a rgid_arr <<< "~{sep=' ' fastqs_se_readgroup_id}"

        # Validate lengths
        len=${#fq_arr[@]}
        if [[ ${#rg_arr[@]} -ne $len || ${#rgid_arr[@]} -ne $len ]]; then
            echo "ERROR: input array lengths differ" 1>&2
            exit 1
        fi

        # Write header
        printf "fastq\treadgroup\treadgroup_id\n" > records.tsv

        # Emit rows
        for ((i=0; i<len; i++)); do
            printf "%s\t%s\t%s\t%s\n" \
                "${fq_arr[i]}" \
                "${rg_arr[i]}" \
                "${rgid_arr[i]}" >> records.tsv
        done
    >>>

    output {
        Array[FastqSingleRecord] fastq_single_records = read_objects("records.tsv")
    }

    runtime {
        docker: "~{docker_image}~{docker_image_hash_or_tag}"
        memory: "~{mem_gb} GB"
        disks: "local-disk ~{disk_space} SSD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }
}

task bwa_pe {
    input {
        FastqPairRecord fastq_record
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        File ref_amb
        File ref_ann
        File ref_bwt
        File ref_pac
        File ref_sa
        File? ref_alt
        Int cpu = 8
        Int preemptible = 2
        Int max_retries = 1
        Int additional_memory_mb = 0
        Int additional_disk_gb = 0
    }

    File fastq1 = fastq_record.forward_fastq
    File fastq2 = fastq_record.reverse_fastq
    String readgroup = fastq_record.readgroup
    String outbam = fastq_record.readgroup_id + ".bam"
    Float ref_size = size([ref_fasta, ref_dict, ref_amb, ref_ann, ref_bwt, ref_pac, ref_sa, ref_fasta_index], "GiB")
    Int computed_mem = 8000 + ceil(size([fastq1, fastq2], "MiB") * 0.1) + additional_memory_mb
    Int mem = if computed_mem < 16000 then computed_mem else 16000
    Int disk_space = ceil((size([fastq1, fastq2], "GiB") * 2) + ref_size) + 10 + additional_disk_gb

    command <<<
        set -euo pipefail

        bwa mem -K 100000000 \
            -t ~{cpu} \
            -T 0 \
            -R "~{readgroup}" \
            ~{ref_fasta} \
            ~{fastq1} \
            ~{fastq2} \
        | samtools view \
            -@ ~{cpu} \
            -Shb \
            -o ~{outbam} \
            -
    >>>

    output {
        File bam = "~{outbam}"
    }

    runtime {
        docker: "broadgdac/bwa:0.7.15-r1142-dirty"
        memory: mem + " MiB"
        disks: "local-disk " + disk_space + " SSD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }
}

task bwa_se {
    input {
        FastqSingleRecord fastq_record
        File ref_fasta
        File ref_dict
        File ref_amb
        File ref_ann
        File ref_bwt
        File ref_pac
        File ref_sa
        File ref_fasta_index
        File? ref_alt
        Int cpu = 8
        Int preemptible = 2
        Int max_retries = 1
        Int additional_memory_mb = 0
        Int additional_disk_gb = 0
    }

    File fastq = fastq_record.fastq
    String readgroup = fastq_record.readgroup
    String outbam = fastq_record.readgroup_id + ".bam"
    Float ref_size = size([ref_fasta, ref_dict, ref_amb, ref_ann, ref_bwt, ref_pac, ref_sa, ref_fasta_index], "GiB")
    Int computed_mem = 8000 + ceil(size(fastq, "MiB") * 0.1) + additional_memory_mb
    Int mem = if computed_mem < 16000 then computed_mem else 16000
    Int disk_space = ceil((size(fastq, "GiB") * 2) + ref_size) + 10 + additional_disk_gb

    command <<<
        set -euo pipefail

        bwa mem -K 100000000 \
            -t ~{cpu} \
            -T 0 \
            -R "~{readgroup}" \
            ~{ref_fasta} \
            ~{fastq} \
        | samtools view \
            -@ ~{cpu} \
            -Shb \
            -o ~{outbam} \
            -
    >>>

    output {
        File bam = "~{outbam}"
    }

    runtime {
        docker: "broadgdac/bwa:0.7.15-r1142-dirty"
        memory: mem + " MiB"
        disks: "local-disk " + disk_space + " SSD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }
}

task picard_markduplicates {
    input {
        Array[File]+ bams
        String outbam

        Int compression_level = 2
        Int preemptible_tries = 1
        Int max_retries = 1
        String validation_stringency = "SILENT"
        String assume_sort_order = "queryname"

        # The program default for READ_NAME_REGEX is appropriate in nearly every case.
        # Sometimes we wish to supply "null" in order to turn off optical duplicate detection
        # This can be desirable if you don't mind the estimated library size being wrong and optical duplicate detection is taking >7 days and failing
        String? read_name_regex
        Int memory_multiplier = 1
        Int additional_disk_gb = 20

        Float? sorting_collection_size_ratio
    }

    Float total_input_size = size(bams, "GiB")
    String metrics_filename = outbam + ".metrics"
    # The merged bam will be smaller than the sum of the parts so we need to account for the unmerged inputs and the merged output.
    # Mark Duplicates takes in as input readgroup bams and outputs a slightly smaller aggregated bam. Giving .25 as wiggleroom
    Float md_disk_multiplier = 3
    Int disk_size = ceil(md_disk_multiplier * total_input_size) + additional_disk_gb

    Float memory_size = 7.5 * memory_multiplier
    Int java_memory_size = (ceil(memory_size) - 2)

    # Task is assuming query-sorted input so that the Secondary and Supplementary reads get marked correctly
    # This works because the output of BWA is query-grouped and therefore, so is the output of MergeBamAlignment.
    # While query-grouped isn't actually query-sorted, it's good enough for MarkDuplicates with ASSUME_SORT_ORDER="queryname"

    command <<<
        set -euo pipefail

        java -Dsamjdk.compression_level=~{compression_level} -Xms~{java_memory_size}g -jar /usr/picard/picard.jar \
            MarkDuplicates \
            INPUT=~{sep=' INPUT=' bams} \
            OUTPUT=~{outbam} \
            METRICS_FILE=~{metrics_filename} \
            VALIDATION_STRINGENCY=~{validation_stringency} \
            ASSUME_SORT_ORDER=~{assume_sort_order} \
            ~{"SORTING_COLLECTION_SIZE_RATIO=" + sorting_collection_size_ratio} \
            ~{"READ_NAME_REGEX=" + read_name_regex}
    >>>

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
        preemptible: preemptible_tries
        memory: "~{memory_size} GiB"
        disks: "local-disk " + disk_size + " SSD"
    }

    output {
        File bam = "~{outbam}"
        File metrics = "~{metrics_filename}"
    }
}

task sort_and_index_markdup_bam {
    input {
        File input_bam
        String output_bam_basename
        String tmp_prefix = "tmp_srt"
        Int cpu = 8
        Int preemptible = 2
        Int max_retries = 1
        Int additional_memory_mb = 0
        Int additional_disk_gb = 0
    }

    Int computed_mem = ceil(size(input_bam, "MiB") * 0.25) + 12000 + additional_memory_mb
    Int mem = if computed_mem < 48000 then computed_mem else 48000
    Int disk_space = ceil(size(input_bam, "GiB") * 3) + 10 + additional_disk_gb
    Int mem_per_thread = floor(mem / cpu * 0.85)
    Int index_threads = cpu
    String output_bam = output_bam_basename + ".bam"
    String output_bai = output_bam_basename + ".bai"

    command <<<
        set -euo pipefail

        samtools sort \
            -@ ~{cpu} \
            -o ~{output_bam} \
            -T ~{tmp_prefix} \
            -m ~{mem_per_thread}M \
            ~{input_bam}

        samtools index \
            -@ ~{index_threads} \
            -b \
            ~{output_bam} \
            ~{output_bai}
    >>>

    output {
        File output_bam = "~{output_bam}"
        File output_bai = "~{output_bai}"
    }

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/samtools:1.10"
        memory: mem + " MiB"
        disks: "local-disk " + disk_space + " SSD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }
}

task CalculateSomaticContamination {
    input {
        File? intervals
        File reference
        File reference_dict
        File reference_index
        File tumor_cram_or_bam
        File tumor_crai_or_bai
        File? normal_cram_or_bam
        File? normal_crai_or_bai
        File? contamination_vcf
        File? contamination_vcf_index

        # runtime
        String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.2.4.1"
        File? gatk_override
        Int? additional_disk
        Int memory_mb = 3000
        Int? preemptible_attempts
        Int? max_retries
    }

    Int disk_size = ceil(size(tumor_cram_or_bam,"GB") + size(normal_cram_or_bam,"GB")) + select_first([additional_disk, 10])
    Int command_mem = memory_mb - 500

    command <<<
        set -euo pipefail

        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

        java -Xmx~{command_mem}m -jar ${GATK_LOCAL_JAR} GetPileupSummaries \
             -R ~{reference} \
             -I ~{tumor_cram_or_bam} \
             ~{"--interval-set-rule INTERSECTION -L " + intervals} \
             -V ~{contamination_vcf} \
             -L ~{contamination_vcf} \
             -O pileups.table

        if [[ -f "~{normal_cram_or_bam}" ]];
        then
            java -Xmx~{command_mem}m -jar ${GATK_LOCAL_JAR} GetPileupSummaries \
                 -R ~{reference} \
                 -I ~{normal_cram_or_bam} \
                 ~{"--interval-set-rule INTERSECTION -L " + intervals} \
                 -V ~{contamination_vcf} \
                 -L ~{contamination_vcf} \
                 -O normal_pileups.table

            java -Xmx~{command_mem}m -jar ${GATK_LOCAL_JAR} CalculateContamination \
                 -I pileups.table \
                 -O contamination.table \
                 --tumor-segmentation segments.table \
                 -matched normal_pileups.table
        else
            touch normal_pileups.table
            java -Xmx~{command_mem}m -jar ${GATK_LOCAL_JAR} CalculateContamination \
                 -I pileups.table \
                 -O contamination.table \
                 --tumor-segmentation segments.table
        fi

        grep -v ^sample contamination.table | awk '{print($2)}' > contam.txt
    >>>

    runtime {
        docker: gatk_docker
        bootDiskSizeGb: 12
        memory: memory_mb + " MiB"
        maxRetries: select_first([max_retries, 2])
        disks: "local-disk " + disk_size + " SSD"
        preemptible: select_first([preemptible_attempts, 3])
    }

    output {
        File pileups = "pileups.table"
        File normal_pileups = "normal_pileups.table"
        File contamination_table = "contamination.table"
        File maf_segments = "segments.table"
        File contamination = "contam.txt"
    }
}

task check_bqsr_and_oq_tags {
    input {
        File bam
        File bai
        File ref_fasta

        String docker_image = "us-central1-docker.pkg.dev/depmap-omics/terra-images/samtools"
        String docker_image_hash_or_tag = ":production"
        Int mem_gb = 2
        Int cpu = 1
        Int preemptible = 3
        Int max_retries = 0
        Int additional_disk_gb = 0
    }

    Int disk_space = ceil(size(bam, "GiB")) + 10 + additional_disk_gb

    parameter_meta {
        bam: { localization_optional: true }
        bai: { localization_optional: true }
        ref_fasta: { localization_optional: true }
    }

    command <<<
        set -euo pipefail

        # set up auth for samtools streaming from GCS buckets
        export GCS_OAUTH_TOKEN="$(gcloud auth application-default print-access-token)"
        export GCS_REQUESTER_PAYS_PROJECT="$(gcloud config get-value project -q)"

        # get header (and single read) to check for previous BQSR recalibration
        samtools head "~{bam}" -n 1 --reference "~{ref_fasta}" > head.txt

        BQSR_PERFORMED=$(grep -q "ApplyBQSR" head.txt && echo "true" || echo "false") && \
            echo "${BQSR_PERFORMED}" > bqsr_performed.txt

        if [ "${BQSR_PERFORMED}" = "true" ]; then
            echo "BQSR has already been performed on ~{bam}"
        else
            echo "BQSR has not been performed on ~{bam} yet"
        fi

        HAS_OQ_TAGS=$(tail -n 1 head.txt | awk '
            BEGIN {
                FS="\t"
                found = 0
            }

            {
                for (i=1; i<=NF; i++) {
                    if ($i ~ /^OQ:/) {
                        found = 1
                        print "true"
                        exit
                    }
                }
            }

            END {
                if (!found) {
                    print "false"
                }
            }
        ') && echo "${HAS_OQ_TAGS}" > has_oq_tags.txt

        if [ "${HAS_OQ_TAGS}" = "true" ]; then
            echo "Found OQ tags; BQSR can be performed again with --use-original-qualities"
        else
            echo "Reads are missing OQ tag; BQSR cannot be performed again"
        fi
    >>>

    output {
        Boolean bqsr_performed = read_boolean("bqsr_performed.txt")
        Boolean has_oq_tags = read_boolean("has_oq_tags.txt")
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

# BQSR scatter-gather tasks from
# https://github.com/broadinstitute/warp/blob/9942eed7d6c8aba87ddf096fc77cddc6f5052e54/tasks/broad/Utilities.wdl
task CreateSequenceGroupingTSV {
    input {
        File ref_dict

        Int mem_gb = 1
        Int cpu = 1
        Int preemptible = 3
        Int max_retries = 0
        Int additional_disk_gb = 0
    }

    Int disk_space = 1 + additional_disk_gb

    # Use python to create the Sequencing Groupings used for BQSR and PrintReads Scatter.
    # It outputs to stdout where it is parsed into a wdl Array[Array[String]]
    # e.g. [["1"], ["2"], ["3", "4"], ["5"], ["6", "7", "8"]]
    command <<<
        set -euo pipefail

        python3 <<CODE
        with open("~{ref_dict}", "r") as ref_dict_file:
            sequence_tuple_list = []
            longest_sequence = 0

            for line in ref_dict_file:
                if line.startswith("@SQ"):
                    line_split = line.split("\t")
                    # (Sequence_Name, Sequence_Length)
                    sequence_tuple_list.append((line_split[1].split("SN:")[1], int(line_split[2].split("LN:")[1])))

            longest_sequence = sorted(sequence_tuple_list, key=lambda x: x[1], reverse=True)[0][1]

        # We are adding this to the intervals because hg38 has contigs named with embedded colons and a bug in GATK strips off
        # the last element after a :, so we add this as a sacrificial element.
        hg38_protection_tag = ":1+"
        # initialize the tsv string with the first sequence
        tsv_string = sequence_tuple_list[0][0] + hg38_protection_tag
        temp_size = sequence_tuple_list[0][1]

        for sequence_tuple in sequence_tuple_list[1:]:
            if temp_size + sequence_tuple[1] <= longest_sequence:
                temp_size += sequence_tuple[1]
                tsv_string += "\t" + sequence_tuple[0] + hg38_protection_tag
            else:
                tsv_string += "\n" + sequence_tuple[0] + hg38_protection_tag
                temp_size = sequence_tuple[1]

        # add the unmapped sequences as a separate line to ensure that they are recalibrated as well
        with open("sequence_grouping.txt","w") as tsv_file:
            tsv_file.write(tsv_string)
            tsv_file.close()

        tsv_string += '\n' + "unmapped"

        with open("sequence_grouping_with_unmapped.txt","w") as tsv_file_with_unmapped:
            tsv_file_with_unmapped.write(tsv_string)
            tsv_file_with_unmapped.close()
        CODE
    >>>

    output {
        Array[Array[String]] sequence_grouping = read_tsv("sequence_grouping.txt")
        Array[Array[String]] sequence_grouping_with_unmapped = read_tsv("sequence_grouping_with_unmapped.txt")
    }

    runtime {
        docker: "python:3.12.4-slim"
        memory: "~{mem_gb} GB"
        disks: "local-disk ~{disk_space} SSD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }
}

task BaseRecalibrator {
    input {
        Array[String] bqsr_regions
        Int n_bqsr_splits
        File bam
        File bam_index
        File dbsnp_vcf
        File dbsnp_vcf_index
        File ref_dict
        File ref_fasta
        File ref_fasta_index

        Int cpu = 2
        Int preemptible = 2
        Int max_retries = 1
        Int additional_memory_mb = 0
        Int additional_disk_gb = 0
    }

    String output_grp = basename(bam, ".bam") + "_bqsr.grp"

    Float ref_size = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + size(ref_dict, "GiB")
    Float dbsnp_size = size(dbsnp_vcf, "GiB")
    Int disk_space = ceil((size(bam, "GiB") / n_bqsr_splits) + ref_size + dbsnp_size) + 20

    Int mem = ceil(size(bam, "MiB") / n_bqsr_splits) + 6000 + additional_memory_mb
    Int jvm_mem = mem - 1000
    Int max_heap = mem - 500

    parameter_meta {
        bam: { localization_optional: true }
        bam_index: { localization_optional: true }
        dbsnp_vcf: { localization_optional: true }
        dbsnp_vcf_index: { localization_optional: true }
        ref_dict: { localization_optional: true }
        ref_fasta: { localization_optional: true }
        ref_fasta_index: { localization_optional: true }
    }

    command <<<
        set -euo pipefail

        java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal -Xlog:gc=debug:file=gc_log.log -Xms~{jvm_mem}m -Xmx~{max_heap}m \
            -jar /root/gatk.jar BaseRecalibrator \
            --input ~{bam} \
            --use-original-qualities \
            --known-sites ~{dbsnp_vcf} \
            --reference ~{ref_fasta} \
            --tmp-dir . \
            --output ~{output_grp} \
            --intervals ~{sep=" --intervals " bqsr_regions}
    >>>

    output {
        File bqsr_recal_file = "~{output_grp}"
    }

    runtime {
        docker: "us-central1-docker.pkg.dev/depmap-omics/terra-images/gatk:4.6.1.0"
        memory: mem + " MiB"
        disks: "local-disk " + disk_space + " SSD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }
}

task GatherBqsrReports {
    input {
        Array[File] input_bqsr_reports
        String output_report_filename

        Int mem_gb = 4
        Int cpu = 1
        Int preemptible = 3
        Int max_retries = 0
        Int additional_disk_gb = 0
    }

    Int disk_space = 3 * ceil(size(input_bqsr_reports, "GiB")) + additional_disk_gb

    command <<<
        set -euo pipefail

        java -Xms3000m -Xmx3000m -jar /root/gatk.jar \
            GatherBQSRReports \
            -I ~{sep=' -I ' input_bqsr_reports} \
            -O ~{output_report_filename}
    >>>

    output {
        File output_bqsr_report = "~{output_report_filename}"
    }

    runtime {
        docker: "us-central1-docker.pkg.dev/depmap-omics/terra-images/gatk:4.6.1.0"
        memory: "~{mem_gb} GB"
        disks: "local-disk ~{disk_space} SSD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }
}

task ApplyBQSR {
    input {
        Array[String] bqsr_regions
        Int n_bqsr_splits
        File bam
        File bam_index
        File bqsr_recal_file
        Boolean emit_original_quals = false

        Int cpu = 2
        Int preemptible = 2
        Int max_retries = 1
        Int additional_memory_mb = 0
        Int additional_disk_gb = 0
    }

    Int mem = ceil(size(bam, "MiB") / n_bqsr_splits) + 4000 + additional_memory_mb
    Int jvm_mem = mem - 1000
    Int max_heap = mem - 500
    Int disk_space = ceil((size(bam, "GiB") * 3) / n_bqsr_splits) + 20 + additional_disk_gb

    String output_bam = basename(bam)

    parameter_meta {
        bam: { localization_optional: true }
        bam_index: { localization_optional: true }
    }

    command <<<
        set -euo pipefail

        java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal -Xlog:gc=debug:file=gc_log.log -Xms~{jvm_mem}m -Xmx~{max_heap}m \
            -jar /root/gatk.jar ApplyBQSR \
            --input ~{bam} \
            --bqsr-recal-file ~{bqsr_recal_file} \
            --emit-original-quals ~{emit_original_quals} \
            --use-original-qualities \
            --add-output-sam-program-record \
            --tmp-dir . \
            --output ~{output_bam} \
            --intervals ~{sep=" --intervals " bqsr_regions}
    >>>

    output {
        File recalibrated_bam = "~{output_bam}"
    }

    runtime {
        docker: "us-central1-docker.pkg.dev/depmap-omics/terra-images/gatk:4.6.1.0"
        memory: mem + " MiB"
        disks: "local-disk " + disk_space + " SSD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }
}

task GatherBamFiles {
    input {
        Array[File] input_bams
        String output_bam_basename

        Int mem_gb = 4
        Int cpu = 1
        Int preemptible = 2
        Int max_retries = 1
        Int additional_memory_mb = 0
        Int additional_disk_gb = 0
    }

    Int java_max_memory_mb = ceil(mem_gb * 1000) - 500
    Int java_inital_memory_mb = (mem_gb * 1000) - 1000
    Int disk_space = 2 * ceil(size(input_bams, "GiB")) + 10 + additional_disk_gb

    command <<<
        set -euo pipefail

        java -Xms~{java_inital_memory_mb}m -Xmx~{java_max_memory_mb}m -jar /usr/picard/picard.jar \
            GatherBamFiles \
            INPUT=~{sep=' INPUT=' input_bams} \
            OUTPUT=~{output_bam_basename}.bam \
            CREATE_INDEX=true
    >>>

    output {
        File output_bam = "~{output_bam_basename}.bam"
        File output_bai = "~{output_bam_basename}.bai"
    }

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
        memory: mem_gb + " GiB"
        disks: "local-disk " + disk_space + " SSD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }
}

# Collect quality metrics from the aggregated bam
task collect_insert_size_metrics {
    input {
        File input_bam
        String output_bam_basename

        Int cpu = 1
        Int preemptible = 2
        Int max_retries = 1
        Int additional_memory_mb = 0
        Int additional_disk_gb = 0
    }

    Int mem_mb = ceil(size(input_bam, "GiB")) + 7000 + additional_memory_mb
    Int jvm_mem = mem_mb - 1000
    Int max_heap = mem_mb - 500
    Int disk_size = ceil(size(input_bam, "GiB")) + 20 + additional_disk_gb

    command <<<
        set -euo pipefail

        java -Xms~{jvm_mem}m -Xmx~{max_heap}m -jar /usr/picard/picard.jar \
            CollectInsertSizeMetrics \
            INPUT=~{input_bam} \
            OUTPUT=~{output_bam_basename}.insert_size_metrics \
            HISTOGRAM_FILE=~{output_bam_basename}.insert_size_histogram.pdf
    >>>

    output {
        File insert_size_histogram_pdf = "~{output_bam_basename}.insert_size_histogram.pdf"
        File insert_size_metrics = "~{output_bam_basename}.insert_size_metrics"
    }

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
        memory: mem_mb + " MiB"
        disks: "local-disk " + disk_size + " SSD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }
}
