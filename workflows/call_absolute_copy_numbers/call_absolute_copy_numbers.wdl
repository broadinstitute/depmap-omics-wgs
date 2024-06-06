version 1.0

struct PureCnSolution {
    Float purity
    Float ploidy
    Float contamination
    Boolean flagged
    Boolean curated
    String comment
    Float cin
    Float cin_ploidy_robust
    Float cin_allele_specific
    Float cin_allele_specific_ploidy_robust
    Boolean is_wgd
}

workflow call_absolute_copy_numbers {
    input {
        String workflow_version = "1.0"
        String workflow_url # populate this with the public URL of this script

        String sample_id
        File called_segments
        File called_mutations
        File intervals
    }

    call run_purecn {
        input:
            sample_id = sample_id,
            called_segments = called_segments,
            called_mutations = called_mutations,
            intervals = intervals
    }

    call purecn_postprocess {
        input:
            selected_solution = run_purecn.selected_solution,
            loh = run_purecn.loh
    }

    output {
        File purecn_solutions_pdf = run_purecn.solutions_pdf
        File purecn_chromosomes_pdf = run_purecn.chromosomes_pdf
        File purecn_rds = run_purecn.rds
        File purecn_dnacopy = run_purecn.dnacopy
        File purecn_variants = run_purecn.variants
        File purecn_loh = run_purecn.loh
        File purecn_genes = run_purecn.genes
        File purecn_segmentation = run_purecn.segmentation
        File purecn_local_optima_pdf = run_purecn.local_optima_pdf
        PureCnSolution purecn_solution = purecn_postprocess.solution
    }
}

task run_purecn {
    input {
        String sample_id
        File called_segments
        File called_mutations
        File intervals

        String genome = "hg38"
        Int max_copy_number = 8
        Float min_purity = 0.90
        Float max_purity = 0.99
        String fun_segmentation = "Hclust"
        Int max_segments = 1000
        Int min_total_counts = 20
        String otherArguments = "--post-optimize --model-homozygous"

        Int preemptible = 2
        Int max_retries = 1
        Int cpu = 1
        Int mem_gb = 6
        Int additional_disk_gb = 0
    }

    Int disk_space = ceil(
        size(called_segments, "GiB") + size(called_mutations, "GiB")
    ) + 10 + additional_disk_gb

    command <<<
        Rscript /opt/PureCN/PureCN.R \
            --out="~{sample_id}" \
            --sampleid="~{sample_id}" \
            --seg-file="~{called_segments}" \
            --vcf="~{called_mutations}" \
            --intervals="~{intervals}" \
            --genome="~{genome}" \
            ~{"--max-purity " + max_purity} \
            ~{"--min-purity " + min_purity} \
            ~{"--max-copy-number " + max_copy_number} \
            ~{"--fun-segmentation " + fun_segmentation} \
            ~{"--max-segments " + max_segments} \
            ~{"--min-total-counts " + min_total_counts} \
            ~{otherArguments}
    >>>

    output {
        File solutions_pdf = "~{sample_id}.pdf"
        File chromosomes_pdf = "~{sample_id}_chromosomes.pdf"
        File rds = "~{sample_id}.rds"
        File dnacopy = "~{sample_id}_dnacopy.seg"
        File variants = "~{sample_id}_variants.csv"
        File loh = "~{sample_id}_loh.csv"
        File genes = "~{sample_id}_genes.csv"
        File segmentation = "~{sample_id}_segmentation.pdf"
        File selected_solution = "~{sample_id}.csv"
        File local_optima_pdf = "~{sample_id}_local_optima.pdf"
    }

    runtime {
        docker: "us-docker.pkg.dev/depmap-omics/public/purecn:production"
        memory: "~{mem_gb} GiB"
        disks: "local-disk ~{disk_space} SSD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }
}

task purecn_postprocess {
    input {
        File selected_solution
        File loh

        Int preemptible = 2
        Int max_retries = 1
        Int cpu = 1
        Int mem_gb = 2
        Int disk_space = 10
    }

    command <<<
        poetry run python -m purecn_postprocess \
            --solution="~{selected_solution}" \
            --loh="~{loh}" \
            --out="out.json"
    >>>

    output {
        PureCnSolution solution = read_json("out.json")
    }

    runtime {
        docker: "us-docker.pkg.dev/depmap-omics/public/purecn_postprocess:production"
        memory: "~{mem_gb} GiB"
        disks: "local-disk ~{disk_space} SSD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }
}
