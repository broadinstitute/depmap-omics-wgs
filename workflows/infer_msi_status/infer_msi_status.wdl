version 1.0

workflow infer_msi_status {
    input {
        String sample_id
        File bam
        File bai
    }

    call run_msisensor2 {
        input:
            sample_id = sample_id,
            bam = bam,
            bai = bai
    }

    output {
        Float msisensor2_score = run_msisensor2.msisensor2_score
        File msisensor2_output = run_msisensor2.msisensor2_output
        File msisensor2_output_dist = run_msisensor2.msisensor2_output_dist
        File msisensor2_output_somatic = run_msisensor2.msisensor2_output_somatic
    }
}

task run_msisensor2 {
    input {
        String sample_id
        File bam
        File bai

        String docker_image = "us-central1-docker.pkg.dev/depmap-omics/terra-images/msisensor2"
        String docker_image_hash_or_tag = ":production"
        Int cpu = 1
        Int mem_gb = 8
        Int preemptible = 2
        Int max_retries = 1
        Int additional_disk_gb = 0
    }

    String bam_path = basename(bam)
    Int disk_space = ceil(size(bam, "GiB")) + 10 + additional_disk_gb

    command <<<
        set -euo pipefail

        # MSISensor needs the BAM in the workdir
        mv ~{bam} .
        mv ~{bai} .

        msisensor2 msi \
            -M /msisensor2/models_hg38 \
            -t ~{bam_path} \
            -o ~{sample_id}.msisensor2.output

        awk '
            BEGIN {
                ORS = "\n";
            }

            NR % 2 == 1 {
                # parse header line
                split($0, h, /[[:space:]]+/);
                chrom = h[1];
                pos = h[2];
                left_flank = h[3];
                repeat = h[4];
                right_flank = h[5];
                next;
            }

            NR % 2 == 0 {
                # strip prefix
                sub(/^T:[[:space:]]*/, "", $0);

                # strip trailing whitespace
                sub(/[[:space:]]+$/, "", $0);

                # normalize internal whitespace to commas
                gsub(/[[:space:]]+/, ",", $0);

                # emit NDJSON
                printf "{\"chrom\":\"%s\",\"pos\":%s,\"left_flank\":\"%s\",\"repeat\":\"%s\",\"right_flank\":\"%s\",\"distribution\":[%s]}\n",
                    chrom, pos, left_flank, repeat, right_flank, $0
            }
        ' "~{sample_id}.msisensor2.output_dis" > "~{sample_id}.msisensor2.output_dist.ndjson"

        head -2 ~{sample_id}.msisensor2.output | \
            tail -1 | \
            cut -f3 > \
            ~{sample_id}.msisensor2.score
    >>>

    output {
        Float msisensor2_score = read_float("~{sample_id}.msisensor2.score")
        File msisensor2_output = "~{sample_id}.msisensor2.output"
        File msisensor2_output_dist = "~{sample_id}.msisensor2.output_dist.ndjson"
        File msisensor2_output_somatic = "~{sample_id}.msisensor2.output_somatic"
    }

    runtime {
        docker: "~{docker_image}~{docker_image_hash_or_tag}"
        memory: "~{mem_gb} GiB"
        disks: "local-disk ~{disk_space} SSD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }

    meta {
        allowNestedInputs: true
    }
}
