version 1.0

# Adapted from https://github.com/broadinstitute/warp/blob/GDCWholeGenomeSomaticSingleSample_v1.3.1/pipelines/broad/reprocessing/cram_to_unmapped_bams/CramToUnmappedBams.wdl
#
# Modifications:
#   - Separate CRAM/CRAI vs. BAM/BAI workflow inputs consolidated into singular
#     `input_cram_bam` and `input_crai_bai` ones, with necessary logic to handle either
#     file type.
#   - Use SSDs instead of HDDs for faster (de)localization and higher compute instance
#     availability.
#   - Run `RevertSam` with `--RESTORE_HARDCLIPS false` to deal with invalid data
#     appearing in the `XQ` header from DRAGEN-produced CRAM files.
#
# Use this file as a base and manually implement changes in the upstream code in order
# to retain these modifications.

# If the output_map file is provided, it is expected to be a tab-separated file containing a list of all the read group ids
# found in the input_cram / input_bam and the desired name of the unmapped bams generated for each.
# If the file is not provided, the output names of the unmapped bams will be the read_group_id<unmapped_bam_suffix>
workflow CramToUnmappedBams {
    input {
        String input_type
        File input_cram_bam
        File input_crai_bai
        File? ref_fasta
        File? ref_fasta_index
        File? output_map
        String unmapped_bam_suffix = ".unmapped.bam"
        Int additional_disk = 20
    }

    if (input_type == "CRAM") {
        Float cram_size = size(input_cram_bam, "GiB")

        call CramToBam {
            input:
                ref_fasta = select_first([ref_fasta]),
                ref_fasta_index = select_first([ref_fasta_index]),
                cram_file = input_cram_bam,
                output_basename = basename(input_cram_bam, ".cram"),
                disk_size = ceil(cram_size * 6) + additional_disk
        }
    }

    File input_file = select_first([CramToBam.output_bam, input_cram_bam])
    Float input_size = size(input_file, "GiB")

    if (!defined(output_map)) {
        call GenerateOutputMap {
            input:
                input_bam = input_file,
                unmapped_bam_suffix = unmapped_bam_suffix,
                disk_size = ceil(input_size) + additional_disk
        }
    }

    call SplitUpOutputMapFile {
        input:
            read_group_map_file = select_first([output_map, GenerateOutputMap.output_map])
    }

    scatter (rg_map_file in SplitUpOutputMapFile.rg_to_ubam_file) {
        call SplitOutUbamByReadGroup {
            input:
                input_bam = input_file,
                rg_to_ubam_file = rg_map_file,
                disk_size = ceil(input_size * 2) + additional_disk
        }

        String unmapped_bam_filename = basename(SplitOutUbamByReadGroup.output_bam)

        call RevertSam {
            input:
                input_bam = SplitOutUbamByReadGroup.output_bam,
                output_bam_filename = unmapped_bam_filename,
                disk_size = ceil(input_size * 3) + additional_disk
        }

        call SortSam {
            input:
                input_bam = RevertSam.output_bam,
                output_bam_filename = unmapped_bam_filename
        }

        Float unmapped_bam_size = size(SortSam.output_bam, "GiB")

        call ValidateSamFile {
            input:
                input_bam = SortSam.output_bam,
                report_filename = unmapped_bam_filename + ".validation_report",
                disk_size = ceil(unmapped_bam_size) + additional_disk
        }
    }

    output {
        Array[File] validation_report = ValidateSamFile.report
        Array[File] unmapped_bams = SortSam.output_bam
    }

    meta {
        allowNestedInputs: true
    }
}

task RevertSam {
    input {
        File input_bam
        String output_bam_filename
        Int disk_size
        Int memory_in_MiB = 3000
    }

    Int java_mem = memory_in_MiB - 1000
    Int max_heap = memory_in_MiB - 500

    command <<<
        set -euo pipefail

        java -Xms~{java_mem}m -Xmx~{max_heap}m -jar /usr/picard/picard.jar \
            RevertSam \
            --INPUT ~{input_bam} \
            --OUTPUT ~{output_bam_filename} \
            --VALIDATION_STRINGENCY LENIENT \
            --ATTRIBUTE_TO_CLEAR FT \
            --ATTRIBUTE_TO_CLEAR CO \
            --ATTRIBUTE_TO_CLEAR PA \
            --ATTRIBUTE_TO_CLEAR OA \
            --ATTRIBUTE_TO_CLEAR XA \
            --RESTORE_HARDCLIPS false \
            --SORT_ORDER coordinate
    >>>

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
        disks: "local-disk " + disk_size + " SSD"
        memory: "~{memory_in_MiB} MiB"
        preemptible: 3
    }

    output {
        File output_bam = output_bam_filename
    }
}

# This task is slower than converting straight from cram to bam (note we stream through sam format in between cram and bam)
# This is currently necessary due to a difference in the way the NM tag is calculated in order to produce a valid bam.
task CramToBam {
    input {
        File ref_fasta
        File ref_fasta_index
        File cram_file
        String output_basename
        Int disk_size
        Int cpu = 4
        Int memory_in_MiB = 7000
    }

    Int index_threads = cpu - 1

    command <<<
        set -euo pipefail

        samtools view -@ ~{cpu} -h -T ~{ref_fasta} ~{cram_file} |
        samtools view -@ ~{cpu} -b -o ~{output_basename}.bam -
        samtools index -@ ~{index_threads} -b ~{output_basename}.bam
        mv ~{output_basename}.bam.bai ~{output_basename}.bai
    >>>

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/samtools:1.0.0-1.11-1624651616"
        cpu: cpu
        memory: "~{memory_in_MiB} MiB"
        disks: "local-disk " + disk_size + " SSD"
        preemptible: 2
    }

    output {
        File output_bam = "~{output_basename}.bam"
        File output_bam_index = "~{output_basename}.bai"
    }
}

task GenerateOutputMap {
    input {
        File input_bam
        String unmapped_bam_suffix
        Int disk_size
        Int memory_in_MiB = 3000
    }

    command <<<
        set -euo pipefail

        samtools view -H ~{input_bam} | grep '^@RG' | cut -f2 | sed s/ID:// > readgroups.txt

        echo -e "#READ_GROUP_ID\tOUTPUT" > output_map.tsv

        for rg in `cat readgroups.txt`; do
            echo -e "$rg\t$rg~{unmapped_bam_suffix}" >> output_map.tsv
        done
    >>>

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/samtools:1.0.0-1.11-1624651616"
        disks: "local-disk " + disk_size + " SSD"
        memory: "~{memory_in_MiB} MiB"
        preemptible: 3
    }

    output {
        File output_map = "output_map.tsv"
    }
}

task SplitUpOutputMapFile {
    input {
        File read_group_map_file
        Int disk_size = 10
        Int memory_in_MiB = 3000
    }

    command <<<
        set -euo pipefail

        mkdir rgtemp
        cd rgtemp

        # splits list of mappings into single files.  One line each.
        grep -v '^#' ~{read_group_map_file} | split -l 1 - rg_to_ubam_
    >>>

    runtime {
        docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4:latest"
        disks: "local-disk " + disk_size + " SSD"
        memory: "~{memory_in_MiB} MiB"
    }

    output {
        Array[File] rg_to_ubam_file = glob("rgtemp/rg_to_ubam_*")
    }
}

task SplitOutUbamByReadGroup {
    input {
        File input_bam
        File rg_to_ubam_file
        Int disk_size
        Int cpu = 2
        Int memory_in_MiB = 3000
    }

    Array[Array[String]] tmp = read_tsv(rg_to_ubam_file)

    command <<<
        set -euo pipefail

        echo "Read Group ~{tmp[0][0]} from ~{input_bam} is being written to ~{tmp[0][1]}"
        samtools view -@ ~{cpu} -b -h -r ~{tmp[0][0]} -o ~{tmp[0][1]} ~{input_bam}
    >>>

    output {
        File output_bam = tmp[0][1]
    }

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/samtools:1.0.0-1.11-1624651616"
        cpu: cpu
        disks: "local-disk " + disk_size + " SSD"
        memory: "~{memory_in_MiB} MiB"
        preemptible: 3
    }
}

task ValidateSamFile {
    input {
        File input_bam
        String report_filename
        Int disk_size
        Int memory_in_MiB = 3000
    }

    Int java_mem = memory_in_MiB - 1000
    Int max_heap = memory_in_MiB - 500

    command <<<
        set -euo pipefail

        java -Xms~{java_mem}m -Xmx~{max_heap}m -jar /usr/picard/picard.jar \
            ValidateSamFile \
            --INPUT ~{input_bam} \
            --OUTPUT ~{report_filename} \
            --MODE VERBOSE \
            --IS_BISULFITE_SEQUENCED false
    >>>

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
        disks: "local-disk " + disk_size + " SSD"
        memory: "~{memory_in_MiB} MiB"
        preemptible: 3
    }

    output {
        File report = "~{report_filename}"
    }
}

task SortSam {
    input {
        File input_bam
        String output_bam_filename
        Int memory_in_MiB = 7000
        Float sort_sam_disk_multiplier = 6
    }
    # SortSam spills to disk a lot more because we are only store 300000 records in RAM now because its faster for our data so it needs
    # more disk space.  Also it spills to disk in an uncompressed format so we need to account for that with a larger multiplier
    Int disk_size = ceil(sort_sam_disk_multiplier * size(input_bam, "GiB")) + 20
    Int java_mem = memory_in_MiB - 1000
    Int max_heap = memory_in_MiB - 500

    command <<<
        set -euo pipefail

        java -Xms~{java_mem}m -Xmx~{max_heap}m -jar /usr/picard/picard.jar \
            SortSam \
            --INPUT ~{input_bam} \
            --OUTPUT ~{output_bam_filename} \
            --SORT_ORDER queryname \
            --MAX_RECORDS_IN_RAM 300000
    >>>

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
        disks: "local-disk " + disk_size + " SSD"
        memory: "~{memory_in_MiB} MiB"
        preemptible: 3
    }

    output {
        File output_bam = output_bam_filename
    }
}
