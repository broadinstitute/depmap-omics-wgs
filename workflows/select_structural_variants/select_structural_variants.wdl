version 1.0

workflow select_structural_variants {
    input {
        String sample_id
        File input_bedpe
        File gene_annotation
        File del_annotation
        File dup_annotation
        File cosmic_fusion_gene_pairs
    }

    call do_select_structural_variants {
        input:
            sample_id = sample_id,
            input_bedpe = input_bedpe,
            gene_annotation = gene_annotation,
            del_annotation = del_annotation,
            dup_annotation = dup_annotation,
            cosmic_fusion_gene_pairs = cosmic_fusion_gene_pairs
    }

    output {
        File expanded_bedpe = do_select_structural_variants.expanded_bedpe
        File expanded_filtered_bedpe = do_select_structural_variants.expanded_filtered_bedpe
    }
}

task do_select_structural_variants {
    input {
        String sample_id
        File input_bedpe
        File gene_annotation
        File del_annotation
        File dup_annotation
        File cosmic_fusion_gene_pairs

        String docker_image
        String docker_image_hash_or_tag
        Int mem_gb = 8
        Int cpu = 1
        Int preemptible = 2
        Int max_retries = 0
        Int additional_disk_gb = 0
    }

    Int disk_space = ceil(
        3 * size(input_bedpe, "GiB")
        + size([gene_annotation, del_annotation, dup_annotation], "GiB")
    ) + 1

    command <<<
        set -euo pipefail

        python -m select_structural_variants \
            --input-bedpe="~{input_bedpe}" \
            --gene-annotation="~{gene_annotation}" \
            --del-annotation="~{del_annotation}" \
            --dup-annotation="~{dup_annotation}" \
            --cosmic-fusion-gene-pairs="~{cosmic_fusion_gene_pairs}" \
            --out="~{sample_id}.somatic_structural_variants.parquet"
    >>>

    output {
        File expanded_bedpe = "~{sample_id}.expanded_reannotated.bedpe"
        File expanded_filtered_bedpe = "~{sample_id}_expanded_reannotated_filtered.bedpe"
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
