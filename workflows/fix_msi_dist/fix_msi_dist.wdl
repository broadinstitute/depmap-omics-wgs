version 1.0

workflow fix_msi_dist {
    input {
        String sample_id
        File msisensor2_output_dis
    }

    call do_fix {
        input:
            sample_id = sample_id,
            msisensor2_output_dis = msisensor2_output_dis
    }

    output {
        File msisensor2_output_dist = do_fix.msisensor2_output_dist
    }
}

task do_fix {
    input {
        String sample_id
        File msisensor2_output_dis

        String docker_image = "us-central1-docker.pkg.dev/depmap-omics/terra-images/bedtools2"
        String docker_image_hash_or_tag = ":production"
        Int cpu = 1
        Int mem_gb = 2
        Int preemptible = 2
        Int max_retries = 0
        Int additional_disk_gb = 0
    }

    Int disk_space = 10

    command <<<
        set -euo pipefail

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
                printf(
                    "{\"chrom\":\"%s\",\"pos\":%s,\"left_flank\":\"%s\",\"repeat\":\"%s\",\"right_flank\":\"%s\",\"distribution\":[%s]}\n",
                    chrom, pos, left_flank, repeat, right_flank, $0
                );
            }
        ' "~{msisensor2_output_dis}" > "~{sample_id}.msisensor2.output_dist.ndjson"
    >>>

    output {
        File msisensor2_output_dist = "~{sample_id}.msisensor2.output_dist.ndjson"
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
