version 1.0

workflow annotate_structural_variants {
    input {
        String workflow_version = "1.0"
        String workflow_source_url # populated automatically with URL of this script

        String sample_id
        File vcf

        # pre/post variant filtering
        String include_string_pre_annot
        String include_string_post_annot

        # scatter/gather VCF by chromosome
        File xy_intervals

        # Ensembl VEP annotation
        File ref_fasta_bgz
        File gnomad
        File gnomad_idx
        String vep_chrom_cache_url_prefix
        String vep_chrom_cache_url_suffix

        # reannotation and final processing
        File gtf_bed
        File cosmic_fusion_gene_pairs
    }

    call filter_variants as filter_variants_pre_annot {
        input:
            sample_id = sample_id,
            vcf = vcf,
            include_string = include_string_pre_annot
    }

    call split_vcf_by_chrom {
        input:
            vcf = filter_variants_pre_annot.vcf_filtered,
            xy_intervals = xy_intervals
    }

    scatter (vcf in split_vcf_by_chrom.vcfs) {
        String chrom_num = sub(sub(basename(vcf), "^chr", ""), ".vcf.gz$", "")
        File vep_cache = vep_chrom_cache_url_prefix + chrom_num + vep_chrom_cache_url_suffix

        call ensembl_vep {
            input:
                vcf = vcf,
                output_file_base_name = sample_id + "_chr" + chrom_num + "_ensembl_vep_annot",
                vep_cache = vep_cache,
                ref_fasta_bgz = ref_fasta_bgz,
                gnomad = gnomad,
                gnomad_idx = gnomad_idx
        }
    }

    call gather_vcfs as ensembl_vep_gathered {
        input:
            vcfs = ensembl_vep.vcf_annot,
            output_file_base_name = sample_id + "_ensembl_vep_annot",
    }

    call filter_variants as filter_variants_post_annot {
        input:
            sample_id = sample_id,
            vcf = ensembl_vep_gathered.output_vcf,
            include_string = include_string_post_annot
    }

    call convert_to_bedpe {
        input:
            sample_id = sample_id,
            vcf = filter_variants_post_annot.vcf_filtered
    }

    call reannotate_genes {
        input:
            input_bedpe = convert_to_bedpe.output_bedpe,
            sample_id = sample_id,
            gtf_bed = gtf_bed
    }

    call process_sv_bedpe {
        input:
            sample_id = sample_id,
            input_bedpe = convert_to_bedpe.output_bedpe,
            gene_annotation = reannotate_genes.output_reannotated_bedpe,
            del_annotation = reannotate_genes.annotated_overlap_del,
            dup_annotation = reannotate_genes.annotated_overlap_dup,
            cosmic_fusion_gene_pairs = cosmic_fusion_gene_pairs
    }

    output {
        File sv_annot_vcf = ensembl_vep_gathered.output_vcf
        File sv_annot_bedpe = convert_to_bedpe.output_bedpe
        File sv_annot_reannotated_bedpe = reannotate_genes.output_reannotated_bedpe
        File sv_annot_expanded_bedpe = process_sv_bedpe.expanded_bedpe
        File sv_annot_expanded_filtered_sv_bedpe = process_sv_bedpe.expanded_filtered_bedpe
    }
}

task filter_variants {
    input {
        String sample_id
        File vcf
        String include_string

        String docker_image
        String docker_image_hash_or_tag
        Int mem_gb = 4
        Int cpu = 1
        Int preemptible = 2
        Int max_retries = 1
        Int additional_disk_gb = 0
    }

    Int disk_space = (
        ceil(2 * size(vcf, "GiB")) + 10 + additional_disk_gb
    )

    command <<<
        set -euo pipefail

        bcftools view \
            "~{vcf}" \
            --include='~{include_string}' \
            --output="~{sample_id}_filtered.vcf.gz"
    >>>

    output {
        File vcf_filtered = "~{sample_id_filtered.vcf.gz"
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
        Int preemptible = 2
        Int max_retries = 1
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

task ensembl_vep {
    input {
        File vcf
        String output_file_base_name
        File vep_cache
        File ref_fasta_bgz
        File gnomad
        File gnomad_idx
        Int max_sv_size = 50000000

        String docker_image
        String docker_image_hash_or_tag
        Int mem_gb = 8
        Int cpu = 2
        Int preemptible = 2
        Int max_retries = 1
        Int additional_disk_gb = 0
    }

    Int disk_space = (
        ceil(10 * size(vcf, "GiB") +
             size(ref_fasta_bgz, "GiB") +
             size(gnomad, "GiB") +
             3 * size(vep_cache, "GiB")
        ) + 10 + additional_disk_gb
    )

    command <<<
        set -euo pipefail

        echo "Populating VEP cache data"
        mkdir -p "vep_cache/homo_sapiens"
        tar -C "vep_cache/homo_sapiens" -xzf "~{vep_cache}"

        echo "Running VEP"
        /opt/vep/src/ensembl-vep/vep \
            --assembly="GRCh38" \
            --biotype \
            --buffer_size=5000 \
            --cache \
            --canonical \
            --ccds \
            --dir_cache="vep_cache" \
            --dir_plugins="/plugins" \
            --everything \
            --fasta="~{ref_fasta_bgz}" \
            --fork=~{cpu} \
            --hgvs \
            --input_file="~{vcf}" \
            --mane \
            --max_sv_size=~{max_sv_size} \
            --no_stats \
            --numbers \
            --offline \
            --output_file="~{output_file_base_name}.vcf" \
            --pick \
            --plugin="StructuralVariantOverlap,file=/tmp/Plugins/~{gnomad},same_type=1,overlap_cutoff=80,reciprocal=0,same_type=1,fields=AC%AF" \
            --plugin="pLI,pLI.tsv" \
            --protein \
            --shift_hgvs=0 \
            --species="homo_sapiens" \
            --symbol \
            --terms="SO" \
            --total_length \
            --vcf

        # older version of bgzip in the image doesn't have `--output` option
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

task gather_vcfs {
    input {
        Array[File] vcfs
        String output_file_base_name

        String docker_image
        String docker_image_hash_or_tag
        Int mem_gb = 16
        Int cpu = 2
        Int preemptible = 2
        Int max_retries = 1
        Int additional_disk_gb = 0
    }

    Int command_mem_mb = 1000 * mem_gb - 500
    Int disk_space = ceil(3 * size(vcfs, "GiB")) + 10 + additional_disk_gb

    parameter_meta {
        vcfs: { localization_optional: true }
    }

    command <<<
        set -euo pipefail

        java -Xmx~{command_mem_mb}m -jar /root/gatk.jar GatherVcfsCloud \
            --input ~{sep=" --input " vcfs} \
            --output "combined.vcf.gz" \
            --gather-type "BLOCK" \
            --disable-contig-ordering-check

        java -Xmx~{command_mem_mb}m -jar /root/gatk.jar SortVcf \
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

task convert_to_bedpe {
    input {
        String sample_id
        File vcf

        Int mem_gb = 4
        Int cpu = 1
        Int preemptible = 2
        Int max_retries = 1
        Int additional_disk_gb = 0
    }

    Int disk_space = ceil(10 * size(vcf, "GiB")) + 10

    command <<<
        set -euo pipefail

        VcfToBedpe -in "~{vcf}" -out "~{sample_id}.bedpe"
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
        File output_bedpe = "~{sample_id}.bedpe"
    }

    meta {
        allowNestedInputs: true
    }
}

task reannotate_genes {
    # since VEP doesn't correctly annotate genes at breakpoints, we have to intersect
    # them ourselves
    input {
        String sample_id
        File input_bedpe
        File gtf_bed

        Int mem_gb = 8
        Int cpu = 1
        Int preemptible = 2
        Int max_retries = 1
        Int additional_disk_gb = 0
    }

    Int disk_space = ceil(3 * size(input_bedpe, "GiB")) + ceil(size(gtf_bed, "GiB")) + 10

    command <<<
        set -euo pipefail

        # annotate genes at breakpoint A
        sed '/^#/d' "~{input_bedpe}" | \
            cut -f1-3,13,16 | \
            bedtools intersect -a stdin -b "~{gtf_bed}" -wao | \
            sed 's/;//g' | sed 's/"//g' | \
            awk -F"\t" '{
                if ($5 != ".") {
                    split($4,arr,":");
                    print $1"\t"$2"\t"$3"\t"$4"\t"arr[1]":"arr[2]":"arr[3]":"arr[4]":"arr[5]":"arr[6]":"arr[7]"\t.\t"$(NF-1)
                }
                else {
                    print $1"\t"$2"\t"$3"\t"$4"\t"$4"\t.\t"$9
                }
            }' | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | \
            bedtools groupby -g 1,2,3,4,5 -c 7 -o distinct | sort -k5,5 \
            > gene_overlaps.A.bed # collapse all genes to one row

        # annotate genes at breakpoint B
        sed '/^#/d' "~{input_bedpe}" | \
            cut -f4,5,6,13,16 | \
            bedtools intersect -a stdin -b "~{gtf_bed}" -wao | \
            sed 's/;//g' | sed 's/"//g' | \
            awk -F"\t" '{
                if ($5 != ".") {
                    split($4,arr,":");
                    print $1"\t"$2"\t"$3"\t"$4"\t"arr[1]":"arr[2]":"arr[3]":"arr[4]":"arr[5]":"arr[6]":"arr[7]"\t.\t"$(NF-1)
                }
                else {
                    print $1"\t"$2"\t"$3"\t"$4"\t"$4"\t.\t"$9
                }
            }' | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | \
            bedtools groupby -g 1,2,3,4,5 -c 7 -o distinct | sort -k5,5 \
            > gene_overlaps.B.bed # collapse all genes to one row

        # join breakpoint A annotations and breakpoint B annotations
        join -1 5 -2 5 gene_overlaps.A.bed gene_overlaps.B.bed | \
            sed 's/ /\t/g' | \
            cut -f2- > "~{sample_id}.SV.gene_overlaps.txt"

        # subset DEL and DUP variants
        awk -F"\t" '{ if ($11 == "DEL") print }' "~{input_bedpe}" > "~{sample_id}.DEL.bedpe"
        awk -F"\t" '{ if ($11 == "DUP") print }' "~{input_bedpe}" > "~{sample_id}.DUP.bedpe"

        # for DEL, not just look at genes at the breakpoints, but also genes that lie between breakpoints
        # so here intersect (start_a, end_b) with GTF
        sed '/^#/d' "~{sample_id}.DEL.bedpe" | \
            cut -f1,2,6,13,16 | \
            bedtools intersect -a stdin -b "~{gtf_bed}" -wao | \
            sed 's/;//g' | sed 's/"//g' | \
            awk -F"\t" '{ \
                if ($5 != ".") {
                    split($4,arr,":");
                    print $1"\t"$2"\t"$3"\t"$4"\t"arr[1]":"arr[2]":"arr[3]":"arr[4]":"arr[5]":"arr[6]":"arr[7]"\t.\t"$(NF-1)
                }
                else {
                    print $1"\t"$2"\t"$3"\t"$4"\t"$4"\t.\t"$9
                }
            }' | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | \
            bedtools groupby -g 1,2,3,4,5 -c 7 -o distinct | sort -k5,5 \
            > "~{sample_id}.DEL.overlap.bed"

        # for DUP, not just look at genes at the breakpoints, but also genes that lie between breakpoints
        # so here intersect (start_a, end_b) with GTF
        sed '/^#/d'  "~{sample_id}.DUP.bedpe" | \
            cut -f1,2,6,13,16 | \
            bedtools intersect -a stdin -b "~{gtf_bed}" -wao | \
            sed 's/;//g' | sed 's/"//g' | \
            awk -F"\t" '{
                if ($5 != ".") {
                    split($4,arr,":");
                    print $1"\t"$2"\t"$3"\t"$4"\t"arr[1]":"arr[2]":"arr[3]":"arr[4]":"arr[5]":"arr[6]":"arr[7]"\t.\t"$(NF-1)
                }
                else {
                    print $1"\t"$2"\t"$3"\t"$4"\t"$4"\t.\t"$9
                }
            }' | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | \
            bedtools groupby -g 1,2,3,4,5 -c 7 -o distinct | sort -k5,5 \
            > "~{sample_id}.DUP.overlap.bed"
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
        File output_reannotated_bedpe = "~{sample_id}.SV.gene_overlaps.txt"
        File annotated_overlap_del = "~{sample_id}.DEL.overlap.bed"
        File annotated_overlap_dup = "~{sample_id}.DUP.overlap.bed"
    }

    meta {
        allowNestedInputs: true
    }
}

task process_sv_bedpe {
    # see python module for details, but in short, it includes:
    # expanding INFO annotation column, incorporating gene re-annotation generated above,
    # gnomad filering, rescuing, keeping only columns of interest
    input {
        String sample_id
        File input_bedpe
        File gene_annotation
        File del_annotation
        File dup_annotation
        File cosmic_fusion_gene_pairs

        Int mem_gb = 8
        Int cpu = 1
        Int preemptible = 2
        Int max_retries = 1
        Int additional_disk_gb = 0
    }

    Int disk_space = ceil(
        3 * size(input_bedpe, "GiB")
        + size([gene_annotation, del_annotation, dup_annotation], "GiB")
    ) + 10

    command <<<
        python -u /app/process_sv_bedpe.py \
            "~{input_bedpe}" \
            "~{gene_annotation}" \
            "~{del_annotation}" \
            "~{dup_annotation}" \
            "~{cosmic_fusion_gene_pairs}" \
            "~{sample_id}"
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
        File expanded_bedpe = "~{sample_id}.svs.expanded.reannotated.bedpe"
        File expanded_filtered_bedpe = "~{sample_id}.svs.expanded.reannotated.filtered.bedpe"
    }

    meta {
        allowNestedInputs: true
    }
}
