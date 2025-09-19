version 1.0

workflow call_cnvs {
    input {
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
            coverage = calc_bin_coverage.coverage,
            chrom_sizes = chrom_sizes,
            gc_wig = gc_wig,
            map_wig = map_wig,
            protein_coding_genes_bed = protein_coding_genes_bed,
            ref_fasta = ref_fasta
    }

    output {
        File bin_coverage = calc_bin_coverage.coverage
        File segments = call_segments.segments
        File read_cov_bin = call_segments.read_cov_bin
        File input_params = call_segments.input_params
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

        String docker_image = "us-central1-docker.pkg.dev/depmap-omics/terra-images/samtools_extra"
        String docker_image_hash_or_tag = ":production"
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

        gzip "~{sample_id}.coverage.tsv"
    >>>

    output {
        File coverage = "~{sample_id}.coverage.tsv.gz"
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
        File coverage
        File chrom_sizes
        File gc_wig
        File map_wig
        File protein_coding_genes_bed
        File ref_fasta

        String docker_image = "us-central1-docker.pkg.dev/depmap-omics/terra-images/call_segments"
        String docker_image_hash_or_tag = ":production"
        Int mem_gb = 4
        Int cpu = 1
        Int preemptible = 2
        Int max_retries = 0
        Int additional_disk_gb = 0
    }

    Int disk_space = ceil(
        size(coverage, "GiB")
        + size(gc_wig, "GiB")
        + size(map_wig, "GiB")
        + size(protein_coding_genes_bed, "GiB")
        + size(ref_fasta, "GiB")
    ) + 10 + additional_disk_gb

    command <<<
        set -euo pipefail

        # symlink reference files to working dir
        ln -vs "~{coverage}" "coverage.tsv"
        ln -vs "~{gc_wig}" "gc.wig"
        ln -vs "~{map_wig}" "map.wig"

        # segment with HMM
        Rscript /app/segment.R \
            "coverage.tsv" \
            "gc.wig" \
            "map.wig" \
            "~{sample_id}" \
            100 \
            0.9

        sort \
            --key="1,1V" --key="2,2n" \
            "~{sample_id}.cn_segments.tsv" \
            > "~{sample_id}.cn_segments.sorted.tsv"

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
        ' "~{sample_id}.cn_segments.sorted.tsv" \
        > "~{sample_id}.cn_segments_for_purecn.tsv"

        # collect segment copy numbers for protein-coding genes
        sed '1d' "~{sample_id}.cn_segments.sorted.tsv" \
            | bedtools intersect -b stdin -a "~{protein_coding_genes_bed}" -wao \
            | awk '{ print $1"\t"$2"\t"$3"\t"$4"\t"$11"\t"$12"\t"$11*$12 }' \
            | bedtools groupby -g 1,2,3,4 -c 6,7,7 -o sum,sum,count \
            | awk '{ if ($5 != 0) { print $0"\t"$6/$5 } else { print $0"\tNA" } }' \
            > "~{sample_id}.cn_gene.tsv"

        # calculate weighted-mean copy numbers for protein-coding gene
        sed '1d' "~{sample_id}.read_cov_bin.tsv" \
            | bedtools intersect -a "~{protein_coding_genes_bed}" -b stdin -wao \
            | awk '{ print $0"\t"$17*$18 }' \
            | sort --key="1,1V" --key="2,2n" --key="3,3n" --key="4,4n" \
            | bedtools groupby -g 1,2,3,4 -c 18,19,19 -o sum,sum,count \
            | awk '{ if ($6 == 0) { print $0"\tNA"} else {print $0"\t"$6/($5+1) }}' \
            > "~{sample_id}.cn_gene_weighted_mean.tsv"

        gzip "~{sample_id}.read_cov_bin.tsv"
    >>>

    output {
        File segments = "~{sample_id}.cn_segments_for_purecn.tsv"
        File read_cov_bin = "~{sample_id}.read_cov_bin.tsv.gz"
        File input_params = "~{sample_id}.input_params.tsv"
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
