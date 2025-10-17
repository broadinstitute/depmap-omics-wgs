version 1.0

workflow short_read_remap {
  input {
    File bam_or_cram
    File bam_or_cram_index
    File bed
    File ref_fasta
    File ref_fasta_amb
    File ref_fasta_ann
    File ref_fasta_bwt
    File ref_fasta_pac
    File ref_fasta_sa
    String out_prefix
    Int threads = 1
    File? old_ref  # Required only for CRAM
  }

  call extract_reads {
    input:
      bam_or_cram = bam_or_cram,
      bam_or_cram_index = bam_or_cram_index,
      bed = bed,
      out_prefix = out_prefix,
      threads = threads,
      old_ref = old_ref
  }

  call remap_reads {
    input:
      fq1 = extract_reads.fq1,
      fq2 = extract_reads.fq2,
      ref_fasta = ref_fasta,
      ref_fasta_amb = ref_fasta_amb,
      ref_fasta_ann = ref_fasta_ann,
      ref_fasta_bwt = ref_fasta_bwt,
      ref_fasta_pac = ref_fasta_pac,
      ref_fasta_sa = ref_fasta_sa,
      read_group = extract_reads.read_group,
      out_prefix = out_prefix,
      threads = threads
  }

  output {
    # Outputs from extraction
    File extracted_bam = extract_reads.extracted_bam
    File extracted_sorted_bam = extract_reads.extracted_sorted_bam
    File nonpairs_rnames = extract_reads.nonpairs_rnames
    File pairs_only_bam = extract_reads.pairs_only_bam
    File original_sorted_bam = extract_reads.original_sorted_bam
    File original_sorted_bai = extract_reads.original_sorted_bai
    File fq1 = extract_reads.fq1
    File fq2 = extract_reads.fq2
    File read_group = extract_reads.read_group

    # Outputs from remapping
    File remapped_bam = remap_reads.remapped_bam
    File remapped_sorted_bam = remap_reads.remapped_sorted_bam
    File remapped_sorted_bai = remap_reads.remapped_sorted_bai
  }
}


# =========================
# Task 1: Extract reads
# =========================
task extract_reads {
  input {
    File bam_or_cram
    File bam_or_cram_index
    File bed
    String out_prefix
    Int threads
    File? old_ref
    Int preemptible = 3
  }

  command <<<
    set -euo pipefail

    echo "Extracting reads from: ~{bam_or_cram}"
    echo "Using BED: ~{bed}"
    echo "Threads: ~{threads}"

    if [[ "~{bam_or_cram}" == *.cram ]]; then
      if [ -z "~{old_ref}" ]; then
        echo "Error: CRAM input requires old reference."
        exit 1
      fi
      samtools view -F 2316 -hb --reference "~{old_ref}" -o "~{out_prefix}_extracted_reads.bam" "~{bam_or_cram}" $(perl -ane '{print "$F[0]:$F[1]-$F[2] "}' "~{bed}")
    else
      samtools view -F 2316 -hb -o "~{out_prefix}_extracted_reads.bam" "~{bam_or_cram}" $(perl -ane '{print "$F[0]:$F[1]-$F[2] "}' "~{bed}")
    fi

    # Sort by read name and filter for proper pairs
    samtools view -hb -F 2316 "~{out_prefix}_extracted_reads.bam" | samtools sort -@ "~{threads}" -n -o "~{out_prefix}_extracted_reads_sorted.bam"
    samtools view "~{out_prefix}_extracted_reads_sorted.bam" | cut -f1 | uniq -c | awk '$1!=2 {print $2}' > "~{out_prefix}_nonpairs_rnames.txt"
    samtools view -Sh "~{out_prefix}_extracted_reads_sorted.bam" | fgrep -vf "~{out_prefix}_nonpairs_rnames.txt" | samtools view -hb -o "~{out_prefix}_pairs_only.bam"

    samtools sort -@ "~{threads}" -n -o "~{out_prefix}_original_sorted_by_read_names.bam" "~{out_prefix}_pairs_only.bam"
    samtools sort -@ "~{threads}" -o "~{out_prefix}_original_sorted.bam" "~{out_prefix}_pairs_only.bam"
    samtools index "~{out_prefix}_original_sorted.bam"

    # Extract read group and fastq
    samtools view -H "~{out_prefix}_extracted_reads.bam" | grep "^@RG" | sed 's/        /\\t/g' > "~{out_prefix}_bug_regions.RG.txt"
    samtools fastq -@ "~{threads}" "~{out_prefix}_original_sorted_by_read_names.bam" \
      -1 "~{out_prefix}_extract_1.fastq" \
      -2 "~{out_prefix}_extract_2.fastq" \
      -0 /dev/null -s /dev/null -n

  >>>

  output {
    File extracted_bam = "${out_prefix}_extracted_reads.bam"
    File extracted_sorted_bam = "${out_prefix}_extracted_reads_sorted.bam"
    File nonpairs_rnames = "${out_prefix}_nonpairs_rnames.txt"
    File pairs_only_bam = "${out_prefix}_pairs_only.bam"
    File original_sorted_bam = "${out_prefix}_original_sorted.bam"
    File original_sorted_bai = "${out_prefix}_original_sorted.bam.bai"
    File fq1 = "${out_prefix}_extract_1.fastq"
    File fq2 = "${out_prefix}_extract_2.fastq"
    File read_group = "${out_prefix}_bug_regions.RG.txt"
  }

  runtime {
    cpu: threads
    memory: "16G"
    disks: "local-disk 10 HDD"
    preemptible: preemptible
    # docker: "biocontainers/samtools:v1.17-4-deb_cv1"
    docker: "broadgdac/samtools:1.10"
  }
}


# =========================
# Task 2: Remap reads
# =========================
task remap_reads {
  input {
    File fq1
    File fq2
    File ref_fasta
    File ref_fasta_amb
    File ref_fasta_ann
    File ref_fasta_bwt
    File ref_fasta_pac
    File ref_fasta_sa
    File read_group
    String out_prefix
    Int threads
    Int preemptible = 3
  }

  command <<<
    set -euo pipefail

    echo "Remapping reads with BWA MEM..."
    bwa mem -t "~{threads}" -R "$(head -n 1 ~{read_group})" "~{ref_fasta}" "~{fq1}" "~{fq2}" \
      | samtools view -hb -o "~{out_prefix}_remapped.bam"

    samtools sort -@ "~{threads}" -o "~{out_prefix}_remapped_sorted.bam" "~{out_prefix}_remapped.bam"
    samtools index "~{out_prefix}_remapped_sorted.bam"

  >>>

  output {
    File remapped_bam = "${out_prefix}_remapped.bam"
    File remapped_sorted_bam = "${out_prefix}_remapped_sorted.bam"
    File remapped_sorted_bai = "${out_prefix}_remapped_sorted.bam.bai"
  }

  runtime {
    cpu: threads
    memory: "16G"
    disks: "local-disk 10 SDD"
    # docker: "biocontainers/samtools:v1.17-4-deb_cv1"
    preemptible: preemptible
    docker: "broadgdac/bwa:0.7.15-r1142-dirty"
  }
}