version 1.0

workflow debug_gatk {
    input {
        String workflow_version = "1.1"
        String workflow_source_url # populated automatically with URL of this script
    }

    call check_disk_space

    output {
        Array[File] stats = check_disk_space.stats
    }
}

task check_disk_space {
    command <<<
        du / > du.txt 2>/dev/null
        du -h / > du_h.txt 2>/dev/null
        du --max-depth 0 / > du_0.txt 2>/dev/null
        du --max-depth 1 / > du_1.txt 2>/dev/null
        du -h --max-depth 0 / > du_h_0.txt 2>/dev/null
        du -h --max-depth 1 / > du_h_1.txt 2>/dev/null
        df / > df.txt 2>/dev/null
    >>>

    output {
        Array[File] stats = glob("*.txt")
    }

    runtime {
        docker: "broadinstitute/gatk:4.6.1.0"
        memory: "4 GiB"
        bootDiskSizeGb: 20
        disks: "local-disk 10 SSD"
        preemptible: 0
        maxRetries: 0
        cpu: 1
    }
}
