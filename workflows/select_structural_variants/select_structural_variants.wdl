version 1.0

workflow select_structural_variants {
    input {
        String sample_id
        File input_bedpe
        File gene_annotation
        File del_annotation
        File dup_annotation
        File cosmic_fusion_gene_pairs
        File onco_tsg
    }

    call do_select_structural_variants {
        input:
            sample_id = sample_id,
            input_bedpe = input_bedpe,
            gene_annotation = gene_annotation,
            del_annotation = del_annotation,
            dup_annotation = dup_annotation,
            cosmic_fusion_gene_pairs = cosmic_fusion_gene_pairs,
            onco_tsg = onco_tsg
    }

    output {
        File selected_somatic_sv = do_select_structural_variants.selected_somatic_sv
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
        File onco_tsg
        Int min_depth = 5
        Float sv_gnomad_cutoff = 0.001
        Int large_sv_size = 1000000000

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
            --onco-tsg="~{onco_tsg}" \
            --min_depth=~{min_depth} \
            --sv_gnomad_cutoff=~{sv_gnomad_cutoff} \
            --large-sv-size=~{large_sv_size} \
            --out="~{sample_id}.selected_somatic_sv.parquet"
    >>>

    output {
        File selected_somatic_sv = "~{sample_id}.selected_somatic_sv.parquet"
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
