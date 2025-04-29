version 1.0

workflow annotate_mutations {
    input {
        String workflow_version = "1.0"
        String workflow_source_url # populated automatically with URL of this script

        File ref_fasta
        File ref_fasta_index
        File ref_dict
        String sample_id
        File input_vcf

        # scatter/gather VCF by chromosome
        File? xy_intervals

        # bcftools-based annotations
        Boolean annot_seg_dups = true
        File? segdup_bed
        File? segdup_bed_index

        Boolean annot_repeat_masker = true
        File? repeatmasker_bed
        File? repeatmasker_bed_index

        Boolean annot_hess_drivers = true
        File? hess_drivers
        File? hess_drivers_index

        Boolean annot_oncokb = true
        File? oncokb_annotation

        Boolean annot_civic = true
        File? civic_annotation
        File? civic_annotation_index

        Boolean annot_cosmic_cmc = true
        File? cosmic_cmc
        File? cosmic_cmc_index

        Boolean annot_oncogenes_tsg = true
        File? oncogenes_tsg
        File? oncogenes_tsg_index

        Boolean annot_hgnc = true
        File? hgnc
        File? hgnc_index

        Boolean annot_gc_prop = true
        Int? gc_window

        # snpEff and SnpSift annotation
        Boolean annot_snpeff = true
        Boolean annot_snpsift = true
        File? clinvar_vcf
        File? clinvar_vcf_index

        # open-cravat annotation
        Boolean annot_open_cravat = true
        String? open_cravat_data_sources_url

        # Ensembl VEP annotation
        Boolean annot_ensembl_vep = true
        String vep_pick_order
        File? ref_fasta_bgz
        File? gene_constraint_scores
        File? loftool_scores
        File? alpha_missense
        File? alpha_missense_index
        String? vep_chrom_cache_url_prefix
        String? vep_chrom_cache_url_suffix

        # Funcotator annotation
        Boolean annot_funcotator = true
        String? reference_version = "hg38"
        Boolean? use_gnomad = true
        Boolean? filter_funcotations = false
        File? interval_list
        String? transcript_selection_mode
        Array[String]? transcript_selection_list
        Array[String]? funcotator_annotation_defaults
        Array[String]? funcotator_annotation_overrides
        String? funcotator_data_sources_url
        String? funcotator_extra_args
    }

    if (annot_seg_dups || annot_repeat_masker || annot_hess_drivers || annot_oncokb ||
        annot_civic || annot_cosmic_cmc || annot_oncogenes_tsg || annot_hgnc ||
        annot_gc_prop) {
        call annot_with_bcftools {
            input:
                vcf = input_vcf,
                output_file_base_name = sample_id + "_bcftools_annot",
                annot_seg_dups = annot_seg_dups,
                annot_repeat_masker = annot_repeat_masker,
                annot_hess_drivers = annot_hess_drivers,
                annot_oncokb = annot_oncokb,
                annot_civic = annot_civic,
                annot_cosmic_cmc = annot_cosmic_cmc,
                annot_oncogenes_tsg = annot_oncogenes_tsg,
                annot_hgnc = annot_hgnc,
                annot_gc_prop = annot_gc_prop,
                segdup_bed = segdup_bed,
                segdup_bed_index = segdup_bed_index,
                repeatmasker_bed = repeatmasker_bed,
                repeatmasker_bed_index = repeatmasker_bed_index,
                hess_drivers = hess_drivers,
                hess_drivers_index = hess_drivers_index,
                oncokb_annotation = oncokb_annotation,
                civic_annotation = civic_annotation,
                civic_annotation_index = civic_annotation_index,
                cosmic_cmc = cosmic_cmc,
                cosmic_cmc_index = cosmic_cmc_index,
                oncogenes_tsg = oncogenes_tsg,
                oncogenes_tsg_index = oncogenes_tsg_index,
                hgnc = hgnc,
                hgnc_index = hgnc_index,
                gc_window = gc_window,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index
        }
    }

    if (annot_snpeff || annot_snpsift) {
        call snpeff_snpsift {
            input:
                vcf = input_vcf,
                annot_snpeff = annot_snpeff,
                annot_snpsift = annot_snpsift,
                output_file_base_name = sample_id + "_snpeff_snpsift_annot",
                clinvar_vcf = clinvar_vcf,
                clinvar_vcf_index = clinvar_vcf_index
        }
    }

    if (annot_open_cravat) {
        call open_cravat {
            input:
                vcf = input_vcf,
                output_file_base_name = sample_id + "_open_cravat_annot",
                open_cravat_data_sources_url = select_first([open_cravat_data_sources_url])
        }
    }

    if (annot_ensembl_vep) {
        call split_vcf_by_chrom {
            input:
                vcf = input_vcf,
                xy_intervals = select_first([xy_intervals])
        }

        scatter (vcf in split_vcf_by_chrom.vcfs) {
            String chrom_num = sub(sub(basename(vcf), "^chr", ""), ".vcf.gz$", "")
            File vep_cache = vep_chrom_cache_url_prefix + chrom_num + vep_chrom_cache_url_suffix

            call ensembl_vep {
                input:
                    vcf = vcf,
                    output_file_base_name = sample_id + "_chr" + chrom_num + "_ensembl_vep_annot",
                    vep_cache = vep_cache,
                    vep_pick_order = vep_pick_order,
                    ref_fasta_bgz = select_first([ref_fasta_bgz]),
                    gene_constraint_scores = select_first([gene_constraint_scores]),
                    loftool_scores = select_first([loftool_scores]),
                    alpha_missense = select_first([alpha_missense]),
                    alpha_missense_index = select_first([alpha_missense_index])
            }
        }

        call gather_vcfs as ensembl_vep_gathered {
            input:
                vcfs = ensembl_vep.vcf_annot,
                output_file_base_name = sample_id + "_ensembl_vep_annot",
        }
    }

    if (annot_funcotator) {
        call funcotator {
            input:
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                ref_dict = ref_dict,
                vcf = input_vcf,
                reference_version = select_first([reference_version]),
                output_file_base_name = sample_id + "_funco_annot",
                output_format = "VCF",
                compress = true,
                use_gnomad = select_first([use_gnomad]),
                filter_funcotations = select_first([filter_funcotations]),
                funcotator_data_sources_url = select_first([funcotator_data_sources_url]),
                interval_list = select_first([interval_list]),
                transcript_selection_mode = transcript_selection_mode,
                transcript_selection_list = transcript_selection_list,
                funcotator_annotation_defaults = funcotator_annotation_defaults,
                funcotator_annotation_overrides = funcotator_annotation_overrides,
                extra_args = funcotator_extra_args
        }
    }

    output {
        File? mut_annot_bcftools_vcf = annot_with_bcftools.vcf_annot
        File? mut_annot_snpeff_snpsift_vcf = snpeff_snpsift.vcf_annot
        File? mut_annot_open_cravat_vcf = open_cravat.vcf_annot
        File? mut_annot_vep_vcf = ensembl_vep_gathered.output_vcf
        File? mut_annot_funcotator_vcf = funcotator.vcf_annot
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
        Int preemptible = 3
        Int max_retries = 0
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

task gather_vcfs {
    input {
        Array[File] vcfs
        String output_file_base_name

        String docker_image
        String docker_image_hash_or_tag
        Int mem_gb = 16
        Int cpu = 2
        Int preemptible = 3
        Int max_retries = 0
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

task annot_with_bcftools {
    input {
        File vcf
        String output_file_base_name
        Boolean annot_seg_dups
        Boolean annot_repeat_masker
        Boolean annot_hess_drivers
        Boolean annot_oncokb
        Boolean annot_civic
        Boolean annot_cosmic_cmc
        Boolean annot_oncogenes_tsg
        Boolean annot_hgnc
        Boolean annot_gc_prop
        String? exclude_string
        File? segdup_bed
        File? segdup_bed_index
        File? repeatmasker_bed
        File? repeatmasker_bed_index
        File? hess_drivers
        File? hess_drivers_index
        File? oncokb_annotation
        File? civic_annotation
        File? civic_annotation_index
        File? cosmic_cmc
        File? cosmic_cmc_index
        File? oncogenes_tsg
        File? oncogenes_tsg_index
        File? hgnc
        File? hgnc_index
        Int? gc_window
        File? ref_fasta
        File? ref_fasta_index

        String docker_image
        String docker_image_hash_or_tag
        Int mem_gb = 4
        Int cpu = 1
        Int preemptible = 3
        Int max_retries = 0
        Int additional_disk_gb = 0
    }

    Int disk_space = (
        ceil(size(vcf, "GiB") + 2 * 10 * size(vcf, "GiB")) + 10 + additional_disk_gb
    )

    command <<<
        set -euo pipefail

        TMP_VCF="~{basename(vcf, '.vcf.gz')}_tmp.vcf.gz"

        if ~{annot_seg_dups}; then
            echo "Annotating segmental duplication regions"
            echo '##INFO=<ID=SEGDUP,Number=0,Type=Flag,Description="If variant is in a segmental duplication region">' \
                > segdup.hdr.vcf
            bcftools annotate \
                "~{vcf}" \
                --annotations="~{segdup_bed}" \
                --columns="CHROM,FROM,TO,SEGDUP" \
                --header-lines="segdup.hdr.vcf" \
                --output="${TMP_VCF}"
            rm "~{vcf}" && mv "${TMP_VCF}" "~{vcf}"
        fi

        if ~{annot_repeat_masker}; then
            echo "Annotating repeat masker regions"
            echo '##INFO=<ID=RM,Number=0,Type=Flag,Description="If variant is in a Repeat Masker region">' \
                > repeatmasker.hdr.vcf
            bcftools annotate \
                "~{vcf}" \
                --annotations="~{repeatmasker_bed}" \
                --columns="CHROM,FROM,TO,RM" \
                --header-lines="repeatmasker.hdr.vcf" \
                --output="${TMP_VCF}"
            rm "~{vcf}" && mv "${TMP_VCF}" "~{vcf}"
        fi

        if ~{annot_hess_drivers}; then
            echo "Annotating Hess drivers"

            echo '##INFO=<ID=HESS,Number=1,Type=String,Description="Hess driver signature">' \
                > hess.hdr.vcf

            bcftools annotate \
                "~{vcf}" \
                --annotations="~{hess_drivers}" \
                --columns="CHROM,POS,REF,ALT,HESS" \
                --header-lines="hess.hdr.vcf" \
                --output="${TMP_VCF}"
            rm "~{vcf}" && mv "${TMP_VCF}" "~{vcf}"
        fi

        if ~{annot_oncokb}; then
            echo "Creating OncoKB indexed TSV"
            ONCOKB_TAB_BASENAME="~{basename(select_first([oncokb_annotation]), '.csv')}"

            # get the columns we care about in proper chr-pos sorted order
            csvcut \
                --columns="Chromosome,Position,Ref,Alt,ProteinChange,Oncogenic,MutationEffect,Hotspot" \
                "~{oncokb_annotation}" \
                | csvformat --out-tabs --skip-header \
                | sort --key="1,1V" --key="2,2n" --key="3,3" --key="4,4" \
                > "${ONCOKB_TAB_BASENAME}_to_fix.tsv"

            # URL encode space characters (violates VCF 4.2 spec) and replace True/False
            # with 1/0 in hotspot column
            awk -F "\t" '{
                gsub(/ /, "%20", $6);
                gsub(/ /, "%20", $7);

                if ($8 == "True") {
                    $8 = 1
                } else if ($8 == "False") {
                    $8 = 0
                }

                print
            }' OFS="\t" "${ONCOKB_TAB_BASENAME}_to_fix.tsv" > "${ONCOKB_TAB_BASENAME}.tsv"
            rm "${ONCOKB_TAB_BASENAME}_to_fix.tsv"

            bgzip "${ONCOKB_TAB_BASENAME}.tsv" --output="${ONCOKB_TAB_BASENAME}.tsv.gz"
            tabix "${ONCOKB_TAB_BASENAME}.tsv.gz" -s1 -b2 -e2

            echo "Annotating with OncoKB"

            echo '##INFO=<ID=ONCOKB_PROT,Number=1,Type=String,Description="OncoKB protein change">' \
                > oncokb.hdr.vcf
            echo '##INFO=<ID=ONCOKB_ONCOGENIC,Number=1,Type=String,Description="OncoKB oncogenic">' \
                >> oncokb.hdr.vcf
            echo '##INFO=<ID=ONCOKB_MUTEFF,Number=1,Type=String,Description="OncoKB mutation effect">' \
                >> oncokb.hdr.vcf
            echo '##INFO=<ID=ONCOKB_HOTSPOT,Number=0,Type=Flag,Description="OncoKB hotspot">' \
                >> oncokb.hdr.vcf

            bcftools annotate \
                "~{vcf}" \
                --annotations="${ONCOKB_TAB_BASENAME}.tsv.gz" \
                --columns="CHROM,POS,REF,ALT,ONCOKB_PROT,ONCOKB_ONCOGENIC,ONCOKB_MUTEFF,ONCOKB_HOTSPOT" \
                --header-lines="oncokb.hdr.vcf" \
                --output="${TMP_VCF}"
            rm "~{vcf}" && mv "${TMP_VCF}" "~{vcf}"
        fi

        if ~{annot_cosmic_cmc}; then
            echo "Annotating COSMIC Cancer Mutation Census tiers"

            echo '##INFO=<ID=CMC_TIER,Number=1,Type=Integer,Description="COSMIC CMC Mutation Significance Tier">' \
                > cosmic_cmc.hdr.vcf

            bcftools index "~{vcf}"

            bcftools annotate \
                "~{vcf}" \
                --annotations="~{cosmic_cmc}" \
                --columns="CMC_TIER" \
                --header-lines="cosmic_cmc.hdr.vcf" \
                --output="${TMP_VCF}"
            rm "~{vcf}" && mv "${TMP_VCF}" "~{vcf}"
        fi

        if ~{annot_civic}; then
            echo "Annotating CIVIC variants"

            echo '##INFO=<ID=CIVIC_ID,Number=1,Type=Integer,Description="CIVIC variant ID">' \
                > civic.hdr.vcf
            echo '##INFO=<ID=CIVIC_SCORE,Number=1,Type=Float,Description="CIVIC evidence score">' \
                >> civic.hdr.vcf
            echo '##INFO=<ID=CIVIC_DESC,Number=1,Type=String,Description="CIVIC variant description">' \
                >> civic.hdr.vcf

            bcftools annotate \
                "~{vcf}" \
                --annotations="~{civic_annotation}" \
                --columns="CHROM,POS,REF,ALT,CIVIC_ID,CIVIC_SCORE,CIVIC_DESC" \
                --header-lines="civic.hdr.vcf" \
                --output="${TMP_VCF}"
            rm "~{vcf}" && mv "${TMP_VCF}" "~{vcf}"
        fi

        if ~{annot_oncogenes_tsg}; then
            echo "Annotating oncogenes and tumor suppressor genes"

            echo '##INFO=<ID=ONCOGENE,Number=0,Type=Flag,Description="Is oncogene">' \
                > oncogenes_tsg.hdr.vcf
            echo '##INFO=<ID=TSG,Number=0,Type=Flag,Description="Is tumor suppressor gene">' \
                >> oncogenes_tsg.hdr.vcf

            bcftools annotate \
                "~{vcf}" \
                --annotations="~{oncogenes_tsg}" \
                --columns=CHROM,BEG,END,ONCOGENE,TSG \
                --merge-logic="ONCOGENE:first,TSG:first" \
                --header-lines="oncogenes_tsg.hdr.vcf" \
                --output="${TMP_VCF}"
            rm "~{vcf}" && mv "${TMP_VCF}" "~{vcf}"
        fi

        if ~{annot_hgnc}; then
            echo "Annotating from HGNC"

            echo '##INFO=<ID=HGNC_NAME,Number=.,Type=String,Description="HGNC approved name">' \
                > hgnc.hdr.vcf
            echo '##INFO=<ID=HGNC_GROUP,Number=.,Type=String,Description="HGNC gene group name">' \
                >> hgnc.hdr.vcf

            bcftools annotate \
                "~{vcf}" \
                --annotations="~{hgnc}" \
                --columns=CHROM,BEG,END,HGNC_NAME,HGNC_GROUP \
                --merge-logic="HGNC_NAME:unique,HGNC_GROUP:unique" \
                --header-lines="hgnc.hdr.vcf" \
                --output="${TMP_VCF}"
            rm "~{vcf}" && mv "${TMP_VCF}" "~{vcf}"
        fi

        if ~{annot_gc_prop}; then
            echo "Annotating GC proportion"

            echo '##INFO=<ID=GC_PROP,Number=1,Type=Float,Description="GC proportion in centered ~{gc_window}bp window">' \
                > gc_prop.hdr.vcf

            /app/gc_prop_in_window \
                --vcf="~{vcf}" \
                --fasta="~{ref_fasta}" \
                --output="gc_prop.tsv" \
                --window=~{gc_window}

            tail -n +2 "gc_prop.tsv" > gc_prop.nohead.tsv
            bgzip "gc_prop.nohead.tsv" --output="gc_prop.nohead.tsv.gz"
            tabix "gc_prop.nohead.tsv.gz" -s1 -b2 -e2

            bcftools annotate \
                "~{vcf}" \
                --annotations="gc_prop.nohead.tsv.gz" \
                --columns="CHROM,POS,REF,ALT,GC_PROP" \
                --header-lines="gc_prop.hdr.vcf" \
                --output="${TMP_VCF}"
            rm "~{vcf}" && mv "${TMP_VCF}" "~{vcf}"
        fi

        # ensure it's bgzipped
        bcftools view \
            "~{vcf}" \
            --output="~{output_file_base_name}.vcf.gz" \
            --no-version
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

task snpeff_snpsift {
    input {
        File vcf
        Boolean annot_snpeff
        Boolean annot_snpsift
        String output_file_base_name
        File? clinvar_vcf
        File? clinvar_vcf_index

        String docker_image
        String docker_image_hash_or_tag
        Int mem_gb = 16
        Int cpu = 1
        Int preemptible = 3
        Int max_retries = 0
        Int additional_disk_gb = 0
    }

    Int disk_space = ceil(2 * 10 * size(vcf, "GiB")) + 10 + additional_disk_gb

    command <<<
        set -euo pipefail

        if ~{annot_snpeff}; then
            echo "Annotating with snpEff"
            java -Xmx14g -jar /app/snpEff.jar \
                ann \
                -noStats \
                GRCh38.mane.1.2.ensembl \
                "~{vcf}" \
                > "snpeff_out.vcf"

            bgzip "snpeff_out.vcf" -o "snpeff_out.vcf.gz"
            rm "snpeff_out.vcf" && mv "snpeff_out.vcf.gz" "~{vcf}"
        fi

        if ~{annot_snpsift}; then
            echo "Annotating with SnpSift"
            java -Xmx14g -jar /app/SnpSift.jar \
                annotate \
                -tabix \
                -noDownload \
                ~{clinvar_vcf} \
                "~{vcf}" \
                > "snpsift_out.vcf"

            bgzip "snpsift_out.vcf" -o "snpsift_out.vcf.gz"
            rm "snpsift_out.vcf" && mv "snpsift_out.vcf.gz" "~{vcf}"
        fi

        mv "~{vcf}" "~{output_file_base_name}.vcf.gz"
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

task open_cravat {
    input {
        File vcf
        String output_file_base_name
        String open_cravat_data_sources_url
        Array[String] annotators_to_use = [
            "brca1_func_assay",
            "ccre_screen",
            "gtex",
            "gwas_catalog",
            "pharmgkb",
            "provean",
            "revel",
            "spliceai"
        ]
        String genome = "hg38"
        String modules_options = "vcfreporter.type=separate"

        String docker_image
        String docker_image_hash_or_tag
        Int mem_gb = 16
        Int cpu = 4
        Int preemptible = 3
        Int max_retries = 0
        Int additional_disk_gb = 20
    }

    Int datasources_size_gb = 25
    Int extracted_vcf_size = ceil(2 * 10 * size(vcf, "GiB"))
    Int boot_disk_space = datasources_size_gb + 10
    Int disk_space = extracted_vcf_size + 10 + additional_disk_gb

    command <<<
        set -euo pipefail

        RSYNC_DEST="$(oc config md)"

        echo "Downloading annotator data"
        gcloud storage rsync \
            --recursive \
            --no-user-output-enabled \
            "~{open_cravat_data_sources_url}" \
            $RSYNC_DEST

        echo "Annotating with Open-Cravat"
        oc run ~{vcf} \
            --liftover ~{genome} \
            --mp ~{cpu} \
            --module-option vcfreporter.type=separate \
            -t vcf \
            -d out \
            -a ~{sep=" " annotators_to_use}

        # fix other headers that are difficult for vcf2maf to parse
        sed -e '/^##INFO=<ID=OC_provean__prediction/s/"D(amaging)"/D(amaging)/' \
            -e '/^##INFO=<ID=OC_provean__prediction/s/"N(eutral)"/N(eutral)/' \
            "out/~{basename(vcf)}.vcf" > "~{output_file_base_name}_to_fix.vcf"
        rm "out/~{basename(vcf)}.vcf"

        # URL encode space characters inserted by pharmgkb (violates VCF 4.2 spec)
        awk 'BEGIN {FS="\t"} NF>=5 {gsub(/ /, "%20")} 1' \
            "~{output_file_base_name}_to_fix.vcf" > "~{output_file_base_name}.vcf"
        rm "~{output_file_base_name}_to_fix.vcf"

        bgzip "~{output_file_base_name}.vcf" --threads=~{cpu}
    >>>

    output {
        File vcf_annot = "~{output_file_base_name}.vcf.gz"
    }

    runtime {
        docker: "~{docker_image}~{docker_image_hash_or_tag}"
        memory: "~{mem_gb} GB"
        disks: "local-disk ~{disk_space} SSD"
        bootDiskSizeGb: boot_disk_space
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
        String vep_pick_order
        File gene_constraint_scores
        File loftool_scores
        File alpha_missense
        File alpha_missense_index

        String docker_image
        String docker_image_hash_or_tag
        Int mem_gb = 8
        Int cpu = 2
        Int preemptible = 3
        Int max_retries = 0
        Int additional_disk_gb = 0
    }

    Int disk_space = (
        ceil(10 * size(vcf, "GiB") +
             size(ref_fasta_bgz, "GiB") +
             size(gene_constraint_scores, "GiB") +
             size(loftool_scores, "GiB") +
             size(alpha_missense, "GiB") +
             3 * size(vep_cache, "GiB")
        ) + 10 + additional_disk_gb
    )

    command <<<
        set -euo pipefail

        gunzip -c "~{gene_constraint_scores}" | \
            awk '{ print $2, $20 }' \
            > "pLI.tsv"

        echo "Populating VEP cache data"
        mkdir -p "vep_cache"
        tar -C "vep_cache" -xzf "~{vep_cache}"

        echo "Running VEP"
        /opt/vep/src/ensembl-vep/vep \
            --assembly="GRCh38" \
            --buffer_size=5000 \
            --cache \
            --dir_cache="vep_cache" \
            --dir_plugins="/plugins" \
            --everything \
            --fasta="~{ref_fasta_bgz}" \
            --fork=~{cpu} \
            --input_file="~{vcf}" \
            --no_stats \
            --offline \
            --output_file="~{output_file_base_name}.vcf" \
            --pick \
            --pick-order "~{vep_pick_order}" \
            --plugin="AlphaMissense,file=~{alpha_missense}" \
            --plugin="LoFtool,~{loftool_scores}" \
            --plugin="pLI,pLI.tsv" \
            --species="homo_sapiens" \
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

# Adapted from https://github.com/broadinstitute/gatk/blob/master/scripts/funcotator_wdl/funcotator.wdl
#
# Modifications:
#   - Index the VCF file first
#   - Localize pre-extracted Funcotator datasource from GCS
#   - Use SSD instead of HDD by default
#
# Use this file as a base and manually implement changes in the upstream code in order
# to retain these modifications.
task funcotator {
    input {
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        String funcotator_data_sources_url

        File vcf

        String reference_version

        String output_file_base_name
        String output_format

        Boolean compress
        Boolean use_gnomad

        String? control_id
        String? case_id
        String? sequencing_center
        String? sequence_source
        String? transcript_selection_mode
        Array[String]? transcript_selection_list
        Array[String]? funcotator_annotation_defaults
        Array[String]? funcotator_annotation_overrides
        Array[String]? funcotator_excluded_fields
        Boolean? filter_funcotations
        File? interval_list

        String? extra_args

        String docker_image
        String docker_image_hash_or_tag
        Int mem_gb = 3
        Int cpu = 1
        Int preemptible = 3
        Int max_retries = 0
        Int additional_disk_gb = 0
    }

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = mem_gb * 1000
    Int command_mem = machine_mem - 1000

    # Calculate disk size:
    Float ref_size_gb = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + size(ref_dict, "GiB")
    Float vcf_size_gb = size(vcf, "GiB")
    Float datasources_size_gb = 25
    Int disk_space = ceil(ref_size_gb + datasources_size_gb + vcf_size_gb) + 20

    # Process input args:
    String output_maf = output_file_base_name + ".maf"
    String output_maf_index = output_maf + ".idx"
    String output_vcf = output_file_base_name + if compress then ".vcf.gz" else ".vcf"
    String output_vcf_index = output_vcf +  if compress then ".tbi" else ".idx"
    String output_file = if output_format == "MAF" then output_maf else output_vcf
    String output_file_index = if output_format == "MAF" then output_maf_index else output_vcf_index
    String transcript_selection_arg = if defined(transcript_selection_list) then " --transcript-list " else ""
    String annotation_def_arg = if defined(funcotator_annotation_defaults) then " --annotation-default " else ""
    String annotation_over_arg = if defined(funcotator_annotation_overrides) then " --annotation-override " else ""
    String filter_funcotations_args = if defined(filter_funcotations) && (filter_funcotations) then " --remove-filtered-variants " else ""
    String excluded_fields_args = if defined(funcotator_excluded_fields) then " --exclude-field " else ""
    String interval_list_arg = if defined(interval_list) then " -L " else ""
    String extra_args_arg = select_first([extra_args, ""])

    parameter_meta {
        ref_fasta: { localization_optional: true }
        ref_fasta_index: { localization_optional: true }
        ref_dict: { localization_optional: true }
    }

    command <<<
        set -e
        export GATK_LOCAL_JAR="/root/gatk.jar"
        export GCS_REQUESTER_PAYS_PROJECT="$(gcloud config get-value project -q)"

        # Hack to validate our WDL inputs:
        #
        # NOTE: This happens here so that we don't waste time copying down the data sources if there's an error.
        if [[ "~{output_format}" != "MAF" ]] && [[ "~{output_format}" != "VCF" ]] ; then
            echo "ERROR: Output format must be MAF or VCF."
        fi

        echo "Indexing VCF"
        java -Xmx~{command_mem}m -jar ${GATK_LOCAL_JAR} IndexFeatureFile --input "~{vcf}"

        echo "Downloading Funcotator data sources"
        DATA_SOURCES_FOLDER="./funcotator_data_sources"
        gcloud storage rsync \
            --recursive \
            --no-user-output-enabled \
            ~{funcotator_data_sources_url} \
            $DATA_SOURCES_FOLDER

        if ~{use_gnomad} ; then
            echo "Enabling gnomAD..."
            for potential_gnomad_gz in gnomAD_exome.tar.gz gnomAD_genome.tar.gz ; do
                if [[ -f $DATA_SOURCES_FOLDER/$potential_gnomad_gz ]] ; then
                    cd $DATA_SOURCES_FOLDER
                    tar -zvxf $potential_gnomad_gz
                    cd -
                else
                    echo "ERROR: Cannot find gnomAD folder: $potential_gnomad_gz" 1>&2
                    false
                fi
            done
        fi

        # Run Funcotator:
        java -Xmx~{command_mem}m -jar ${GATK_LOCAL_JAR} Funcotator \
            --data-sources-path $DATA_SOURCES_FOLDER \
            --ref-version ~{reference_version} \
            --output-file-format ~{output_format} \
            -R ~{ref_fasta} \
            -V ~{vcf} \
            -O ~{output_file} \
            --verbosity WARNING \
            ~{interval_list_arg} ~{default="" interval_list} \
            --annotation-default normal_barcode:~{default="Unknown" control_id} \
            --annotation-default tumor_barcode:~{default="Unknown" case_id} \
            --annotation-default Center:~{default="Unknown" sequencing_center} \
            --annotation-default source:~{default="Unknown" sequence_source} \
            ~{"--transcript-selection-mode " + transcript_selection_mode} \
            ~{transcript_selection_arg}~{default="" sep=" --transcript-list " transcript_selection_list} \
            ~{annotation_def_arg}~{default="" sep=" --annotation-default " funcotator_annotation_defaults} \
            ~{annotation_over_arg}~{default="" sep=" --annotation-override " funcotator_annotation_overrides} \
            ~{excluded_fields_args}~{default="" sep=" --exclude-field " funcotator_excluded_fields} \
            ~{filter_funcotations_args} \
            ~{extra_args_arg} \
            --gcs-project-for-requester-pays "${GCS_REQUESTER_PAYS_PROJECT}"

        # Make sure we have a placeholder index for MAF files so this workflow doesn't fail:
        if [[ "~{output_format}" == "MAF" ]] ; then
            touch ~{output_maf_index}
        fi
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
        File vcf_annot = "~{output_file}"
    }
}
