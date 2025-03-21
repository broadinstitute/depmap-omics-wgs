version 1.0

workflow call_cnvs {
    input {
        String workflow_version = "1.0" # internal semver
        String workflow_source_url # populated automatically with URL of this script

        String sample_id
        File bam
        Int min_mapq = 60
        File wgs_bins_bed
        File chrom_sizes
        File gc_wig
        File map_wig
        File protein_coding_genes_bed
        File ref_fasta
    }

    call calc_bin_coverage {
        input:
            sample_id = sample_id,
            bam = bam,
            min_mapq = min_mapq,
            wgs_bins_bed = wgs_bins_bed,
            chrom_sizes = chrom_sizes,
            ref_fasta = ref_fasta
    }

    call call_segments {
        input:
            sample_id = sample_id,
            coverage_wig = calc_bin_coverage.coverage_wig,
            chrom_sizes = chrom_sizes,
            gc_wig = gc_wig,
            map_wig = map_wig,
            protein_coding_genes_bed = protein_coding_genes_bed,
            ref_fasta = ref_fasta
    }

    output {
        File segments = call_segments.segments
        File read_cov_bin = call_segments.read_cov_bin
        File input_params = call_segments.input_params
        File cn_by_gene = call_segments.cn_by_gene
        File cn_by_gene_weighted_mean = call_segments.cn_by_gene_weighted_mean
    }
}

task calc_bin_coverage {
    input {
        String sample_id
        File bam
        Int min_mapq = 60
        File wgs_bins_bed
        File chrom_sizes
        File ref_fasta

        String docker_image
        String docker_image_hash_or_tag
        Int mem_gb = 16
        Int cpu = 4
        Int preemptible = 1
        Int max_retries = 0
        Int additional_disk_gb = 0
    }

    Int disk_space = ceil(
        1.5 * size(bam, "GiB") + size(ref_fasta, "GiB") + size(wgs_bins_bed, "GiB")
    ) + 10 + additional_disk_gb

    command <<<
        set -eu

        echo "collecting raw coverage across provided bins"
        samtools view \
            -T "~{ref_fasta}" \
            -@ ~{cpu} \
            -q ~{min_mapq} \
            --bam "~{bam}" \
        | bedtools coverage \
            -a "~{wgs_bins_bed}" \
            -b stdin \
            -g "~{chrom_sizes}" \
            -sorted \
            -counts \
            -iobuf 300M \
        > "~{sample_id}.coverage.unsorted.bg"

        echo "sorting coverage"
        LC_ALL=C sort \
            -k1,1 \
            -k2,2n \
            "~{sample_id}.coverage.unsorted.bg" \
        > "~{sample_id}.coverage.tsv" && \
        rm "~{sample_id}.coverage.unsorted.bg"

        echo "converting bed graph to bigwig to wig"
        bedGraphToBigWig \
            "~{sample_id}.coverage.tsv" \
            "~{chrom_sizes}" \
            "~{sample_id}.coverage.bw" && \
        rm "~{sample_id}.coverage.tsv"

        bigWigToWig \
            "~{sample_id}.coverage.bw" \
            "~{sample_id}.coverage.wig"
    >>>

    output {
        File coverage_wig = "~{sample_id}.coverage.wig"
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

task call_segments {
    input {
        String sample_id
        File coverage_wig
        File chrom_sizes
        File gc_wig
        File map_wig
        File protein_coding_genes_bed
        File ref_fasta

        String docker_image
        String docker_image_hash_or_tag
        Int mem_gb = 4
        Int cpu = 1
        Int preemptible = 2
        Int max_retries = 0
        Int additional_disk_gb = 0
    }

    Int disk_space = ceil(
        size(coverage_wig, "GiB")
        + size(gc_wig, "GiB")
        + size(map_wig, "GiB")
        + size(protein_coding_genes_bed, "GiB")
        + size(ref_fasta, "GiB")
    ) + 10 + additional_disk_gb

    command <<<
        set -euo pipefail

        # extract coverage column and reformat header lines
        sed "s/#bedGraph section /fixedStep chrom=/g" "~{coverage_wig}" \
            | sed "s/:/ start=/g" \
            | sed "s/\-[0-9]*/ step=1000 span=1000/g" \
            | awk '{ if (index($1, "chr") == 1) { print $4 } else print $0 }' \
            > "~{sample_id}.coverage.formatted.wig"

        # segment with HMM
        Rscript /usr/src/segment.R \
            "~{sample_id}.coverage.formatted.wig" \
            "~{gc_wig}" \
            "~{map_wig}" \
            "~{sample_id}" \
            100 \
            0.9

        # tweak segments file for use as PureCN input
        awk '
            BEGIN {
                FS = OFS = "\t"  # set field separator to tab
            }

            NR == 1 {
                # replace header names
                $1 = "CONTIG"
                $2 = "START"
                $3 = "END"
                $4 = "state"
                $5 = "LOG2_COPY_RATIO_POSTERIOR_50"
                print $0, "NUM_POINTS_COPY_RATIO"
                next
            }

            {
                # calculate NUM_POINTS_COPY_RATIO
                num_points = int(($3 - $2) / 1000)
                print $0, num_points
            }
        ' "~{sample_id}.cn_segments.tsv" \
        > "~{sample_id}.cn_segments_for_purecn.tsv"

        # collect segment copy numbers for protein-coding genes
        sed '1d' "~{sample_id}.cn_segments.tsv" \
            | bedtools intersect -b stdin -a "~{protein_coding_genes_bed}" -wao \
            | awk '{ print $1"\t"$2"\t"$3"\t"$4"\t"$11"\t"$12"\t"$11*$12 }' \
            | bedtools groupby -g 1,2,3,4 -c 6,7,7 -o sum,sum,count \
            | awk '{ if ($5 != 0) { print $0"\t"$6/$5 } else { print $0"\tNA" } }' \
            > "~{sample_id}.cn_gene.tsv"

        # calculate weighted-mean copy numbers for protein-coding gene
        awk '{ if (NR > 1 && $6 >= 0.9) print }' "~{sample_id}.read_cov_bin.tsv" \
            | bedtools intersect -a "~{protein_coding_genes_bed}" -b stdin -wao \
            | awk '{ print $0"\t"$17*$18 }' \
            | sort -k1,1 -k2,2n -k3,3n -k4,4 \
            | bedtools groupby -g 1,2,3,4 -c 18,19,19 -o sum,sum,count \
            | awk '{ print $0"\t"$6/($5+1) }' \
            > "~{sample_id}.cn_gene_weighted_mean.tsv"
    >>>

    output {
        File segments = "~{sample_id}.cn_segments_for_purecn.tsv"
        File read_cov_bin = "~{sample_id}.read_cov_bin.tsv"
        File input_params = "~{sample_id}.input_params.tsv"
        File cn_by_gene = "~{sample_id}.cn_gene.tsv"
        File cn_by_gene_weighted_mean = "~{sample_id}.cn_gene_weighted_mean.tsv"
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
