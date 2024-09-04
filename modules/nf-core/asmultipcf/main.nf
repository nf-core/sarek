process ASMULTIPCF {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-c278c7398beb73294d78639a864352abef2931ce:ba3e6d2157eac2d38d22e62ec87675e12adb1010-0':
        'biocontainers/mulled-v2-c278c7398beb73294d78639a864352abef2931ce:ba3e6d2157eac2d38d22e62ec87675e12adb1010-0' }"

    input:
    tuple val(meta), path(tumor_logr_files)
    tuple val(meta), path(tumor_baf_files)
    tuple val(meta), path(normal_logr_file)
    tuple val(meta), path(normal_baf_file)

    output:
    tuple val(meta), path("*_asmultipcf_segments.txt"), emit: asmultipcf_segments
    tuple val(meta), path("*_asmultipcf_purityploidy.txt"), emit: asmultipcf_purityploidy
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript
    library(ASCAT)

    # Concatenate tumor LogR files
    tumor_logr_data <- do.call(cbind, lapply(strsplit("${tumor_logr_files}", " "), function(file) {
        read.table(file, header = TRUE, check.names = FALSE)
    }))
    write.table(tumor_logr_data, file = "combined_tumor_logr.txt", sep = "\t", quote = FALSE, row.names = FALSE)

    # Concatenate tumor BAF files
    tumor_baf_data <- do.call(cbind, lapply(strsplit("${tumor_baf_files}", " "), function(file) {
        read.table(file, header = TRUE, check.names = FALSE)
    }))
    write.table(tumor_baf_data, file = "combined_tumor_baf.txt", sep = "\t", quote = FALSE, row.names = FALSE)

    # Load the data
    ascat.bc <- ascat.loadData(
        Tumor_LogR_file = "combined_tumor_logr.txt",
        Tumor_BAF_file = "combined_tumor_baf.txt",
        Germline_LogR_file = "$normal_logr_file",
        Germline_BAF_file = "$normal_baf_file"
    )

    # Run multi-sample segmentation
    ascat.bc <- ascat.asmultipcf(ascat.bc, penalty = ${params.ascat_asmultipcf_penalty ?: 5})

    # Run ASCAT
    ascat.output <- ascat.runAscat(ascat.bc)

    # Write out segmented regions
    write.table(ascat.output[["segments"]], file="${prefix}_asmultipcf_segments.txt", sep="\t", quote=FALSE, row.names=FALSE)

    # Write out purity and ploidy info
    purity_ploidy <- data.frame(
        Sample = names(ascat.output\$aberrantcellfraction),
        Purity = unlist(ascat.output\$aberrantcellfraction),
        Ploidy = unlist(ascat.output\$ploidy)
    )
    write.table(purity_ploidy, file="${prefix}_asmultipcf_purityploidy.txt", sep="\t", quote=FALSE, row.names=FALSE)

    # Version export
    writeLines(c("\\"${task.process}\\":", paste0("    ascat: ", packageVersion("ASCAT"))), "versions.yml")
    """
}