version 1.0

workflow compute_molecular_signatures {
    meta {
        description: "Compute molecular signatures using SignatureAnalyzer"
    }

    parameter_meta {
        # inputs
        sample_id: "ID of this sample"
        mut_sig_variants: "Parquet file of variants"
        ref_2bit: "2bit file for the reference genome"
        sa_ref: "SignatureAnalyzer reference signature to use"

        # outputs
        mut_sig_scores: "Parquet file containing the signatures' scores"
    }

    input {
        # per-sample inputs
        String sample_id
        File mut_sig_variants

        # references
        File ref_2bit
        String sa_ref
    }

    call make_maf {
        input:
            sample_id = sample_id,
            mut_sig_variants = mut_sig_variants
    }

    call compute_scores {
        input:
            sample_id = sample_id,
            maf = make_maf.maf,
            ref_2bit = ref_2bit,
            sa_ref = sa_ref
    }

    output {
        File mut_sig_scores = compute_scores.mut_sig_scores
    }
}

task make_maf {
    meta {
        description: "Make a MAF-like file to pass to SignatureAnalyzer"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        sample_id: "ID of this sample"
        mut_sig_variants: "Parquet file of variants"

        # outputs
        maf: "MAF file of the variants"
    }

    input {
        String sample_id
        File mut_sig_variants

        String docker_image  = "us-central1-docker.pkg.dev/depmap-omics/terra-images/compute-molecular-signatures"
        String docker_image_hash_or_tag = ":production"
        Int cpu = 1
        Int mem_gb = 4
        Int preemptible = 2
        Int max_retries = 0
        Int additional_disk_gb = 0
    }

    Int disk_space = ceil(size(mut_sig_variants, "GiB")) * 2 + 10 + additional_disk_gb

    command <<<
        set -euo pipefail

        python -m compute_molecular_signatures \
            make-maf \
            --muts="~{mut_sig_variants}" \
            --maf-out="~{sample_id}.mut_sig_variants.maf.parquet"
    >>>

    output {
        File maf = "~{sample_id}.mut_sig_variants.maf.parquet"
    }

    runtime {
        docker: "~{docker_image}~{docker_image_hash_or_tag}"
        memory: mem_gb + " GiB"
        disks: "local-disk " + disk_space + " SSD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }
}

task compute_scores {
    meta {
        description: "Run SignatureAnalyzer to get the signature scores"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        sample_id: "ID of this sample"
        maf: "MAF file of the variants"
        ref_2bit: "2bit file for the reference genome"
        sa_ref: "SignatureAnalyzer reference signature to use"

        # outputs
        mut_sig_scores: "Parquet file containing the signatures' scores"
    }

    input {
        String sample_id
        File maf
        File ref_2bit
        String sa_ref

        String docker_image  = "us-central1-docker.pkg.dev/depmap-omics/terra-images/signature-analyzer"
        String docker_image_hash_or_tag = ":production"
        Int cpu = 1
        Int mem_gb = 8
        Int preemptible = 1
        Int max_retries = 0
        Int additional_disk_gb = 0
    }

    Int disk_space = ceil(size(maf, "GiB")) + 10 + additional_disk_gb

    command <<<
        set -euo pipefail

        python /app/run_signature_analyzer.py \
            --input-maf="~{maf}" \
            --ref-2bit="~{ref_2bit}" \
            --sa-ref="~{sa_ref}" \
            --parquet-out="~{sample_id}.mut_sig_scores.parquet"
    >>>

    output {
        File mut_sig_scores = "~{sample_id}.mut_sig_scores.parquet"
    }

    runtime {
        docker: "~{docker_image}~{docker_image_hash_or_tag}"
        memory: mem_gb + " GiB"
        disks: "local-disk " + disk_space + " SSD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }
}
