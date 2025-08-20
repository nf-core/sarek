process ASCAT {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4c/4cf02c7911ee5e974ce7db978810770efbd8d872ff5ab3462d2a11bcf022fab5/data'
        : 'community.wave.seqera.io/library/ascat_cancerit-allelecount:c3e8749fa4af0e99'}"

    input:
    tuple val(meta), path(input_normal), path(index_normal), path(input_tumor), path(index_tumor)
    path allele_files
    path loci_files
    path bed_file
    path fasta
    path gc_file
    path rt_file

    output:
    tuple val(meta), path("*alleleFrequencies_chr*.txt"), emit: allelefreqs
    tuple val(meta), path("*BAF.txt"),                    emit: bafs
    tuple val(meta), path("*cnvs.txt"),                   emit: cnvs
    tuple val(meta), path("*LogR.txt"),                   emit: logrs
    tuple val(meta), path("*metrics.txt"),                emit: metrics
    tuple val(meta), path("*png"),                        emit: png
    tuple val(meta), path("*purityploidy.txt"),           emit: purityploidy
    tuple val(meta), path("*segments.txt"),               emit: segments
    path "versions.yml",                                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def gender        = args.gender        ? "${args.gender}"        : "NULL"
    def genomeVersion = args.genomeVersion ? "${args.genomeVersion}" : "NULL"
    def purity        = args.purity        ? "${args.purity}"        : "NULL"
    def ploidy        = args.ploidy        ? "${args.ploidy}"        : "NULL"
    def gc_input      = gc_file            ? "${gc_file}"            : "NULL"
    def rt_input      = rt_file            ? "${rt_file}"            : "NULL"

    def minCounts_arg                   = args.minCounts                   ? ", minCounts = ${args.minCounts}"                                     : ""
    def bed_file_arg                    = bed_file                         ? ", BED_file = '${bed_file}'"                                          : ""
    def chrom_names_arg                 = args.chrom_names                 ? ", chrom_names = ${args.chrom_names}"                                 : ""
    def min_base_qual_arg               = args.min_base_qual               ? ", min_base_qual = ${args.min_base_qual}"                             : ""
    def min_map_qual_arg                = args.min_map_qual                ? ", min_map_qual = ${args.min_map_qual}"                               : ""
    def skip_allele_counting_tumour_arg = args.skip_allele_counting_tumour ? ", skip_allele_counting_tumour = ${args.skip_allele_counting_tumour}" : ""
    def skip_allele_counting_normal_arg = args.skip_allele_counting_normal ? ", skip_allele_counting_normal = ${args.skip_allele_counting_normal}" : ""

    if (args.additional_allelecounter_flags && fasta) {
        additional_allelecounter_arg = ", additional_allelecounter_flags = \"${args.additional_allelecounter_flags} -r ${fasta}\" "
    }
    else if (args.additional_allelecounter_flags) {
        additional_allelecounter_arg = ", additional_allelecounter_flags = \"${args.additional_allelecounter_flags}\" "
    }
    else if (fasta) {
        additional_allelecounter_arg = ", additional_allelecounter_flags = '-r \"${fasta}\"'"
    }
    else {
        additional_allelecounter_arg = ""
    }

    """
    #!/usr/bin/env Rscript
    library(RColorBrewer)
    library(ASCAT)
    options(bitmapType='cairo')

    if(dir.exists("${allele_files}")) {
        # expected production use of a directory
        allele_path   = normalizePath("${allele_files}")
        allele_prefix = paste0(allele_path, "/", "${allele_files}", "_chr")
    } else if(file.exists("${allele_files}")) {
        # expected testing use of a single file
        allele_path   = basename(normalizePath("${allele_files}"))
        allele_prefix = sub('_chr[0-9]+\\\\.txt\$', "_chr", allele_path)
    } else {
        stop("The specified allele files do not exist.")
    }

    if(length(Sys.glob(paste0(allele_prefix,"*")) ) == 0) {
        stop(paste("No allele files found matching", allele_prefix))
    }

    if(dir.exists("${loci_files}")) {
        # expected production use of a directory
        loci_path   = normalizePath("${loci_files}")
        loci_prefix = paste0(loci_path, "/", "${loci_files}", "_chr")
    } else if(file.exists("${loci_files}")) {
        # expected testing use of a single file
        loci_path   = basename(normalizePath("${loci_files}"))
        loci_prefix = sub('_chr[0-9]+\\\\.txt\$', "_chr", loci_path)
    } else {
        stop("The specified loci files do not exist.")
    }

    if(length(Sys.glob(paste0(loci_prefix,"*")) ) == 0) {
        stop(paste("No loci files found matching", loci_prefix))
    }

    # Prepare from BAM files
    ascat.prepareHTS(
        tumourseqfile = "${input_tumor}",
        normalseqfile = "${input_normal}",
        tumourname = paste0("${prefix}", ".tumour"),
        normalname = paste0("${prefix}", ".normal"),
        allelecounter_exe = "alleleCounter",
        alleles.prefix = allele_prefix,
        loci.prefix = loci_prefix,
        gender = "${gender}",
        genomeVersion = "${genomeVersion}",
        nthreads = ${task.cpus}
        ${minCounts_arg}
        ${bed_file_arg}
        ${chrom_names_arg}
        ${min_base_qual_arg}
        ${min_map_qual_arg}
        ${skip_allele_counting_tumour_arg}
        ${skip_allele_counting_normal_arg}
        ${additional_allelecounter_arg}
        , seed = 42
    )

    # Load the data
    ascat.bc = ascat.loadData(
        Tumor_LogR_file = paste0("${prefix}", ".tumour_tumourLogR.txt"),
        Tumor_BAF_file = paste0("${prefix}", ".tumour_tumourBAF.txt"),
        Germline_LogR_file = paste0("${prefix}", ".tumour_normalLogR.txt"),
        Germline_BAF_file = paste0("${prefix}", ".tumour_normalBAF.txt"),
        genomeVersion = "${genomeVersion}",
        gender = "${gender}"
    )

    # Plot the raw data
    ascat.plotRawData(ascat.bc, img.prefix = paste0("${prefix}", ".before_correction."))

    # Optional LogRCorrection
    if("${gc_input}" != "NULL") {

        if(dir.exists("${gc_input}")) {
            # sarek production use of an unzipped folder containing one file
            gc_input = list.files("${gc_input}", recursive = TRUE, full.names = TRUE)
            if(length(gc_input) != 1 | !file.exists(gc_input)) {
                stop("A single gc_input should be provided!")
            }
        } else if(file.exists("${gc_input}")) {
            gc_input = normalizePath("${gc_input}")
        } else {
            stop("gc_input must be a file or folder containing one file")
        }

        if("${rt_input}" != "NULL"){

            if(dir.exists("${rt_input}")) {
                # sarek production use of an unzipped folder containing one file
                rt_input = list.files("${rt_input}", recursive = TRUE, full.names = TRUE)
                if(length(rt_input) != 1 | !file.exists(rt_input)) {
                    stop("A single rt_input should be provided!")
                }
            } else if(file.exists("${rt_input}")) {
                rt_input = normalizePath("${rt_input}")
            } else {
                stop("rt_input must be a file or folder containing one file")
            }

            ascat.bc = ascat.correctLogR(ascat.bc, GCcontentfile = gc_input, replictimingfile = rt_input)
            # Plot raw data after correction
            ascat.plotRawData(ascat.bc, img.prefix = paste0("${prefix}", ".after_correction_gc_rt."))
        }
        else {
            ascat.bc = ascat.correctLogR(ascat.bc, GCcontentfile = gc_input, replictimingfile = ${rt_input})
            # Plot raw data after correction
            ascat.plotRawData(ascat.bc, img.prefix = paste0("${prefix}", ".after_correction_gc."))
        }
    }

    # Segment the data
    ascat.bc = ascat.aspcf(ascat.bc, seed=42)

    # Plot the segmented data
    ascat.plotSegmentedData(ascat.bc)

    # Run ASCAT to fit every tumor to a model, inferring ploidy, normal cell contamination,
    # and discrete copy numbers
    # If psi and rho are manually set:
    if (!is.null(${purity}) && !is.null(${ploidy})){
        ascat.output <- ascat.runAscat(ascat.bc, gamma=1, rho_manual=${purity}, psi_manual=${ploidy})
    } else if(!is.null(${purity}) && is.null(${ploidy})){
        ascat.output <- ascat.runAscat(ascat.bc, gamma=1, rho_manual=${purity})
    } else if(!is.null(${ploidy}) && is.null(${purity})){
        ascat.output <- ascat.runAscat(ascat.bc, gamma=1, psi_manual=${ploidy})
    } else {
        ascat.output <- ascat.runAscat(ascat.bc, gamma=1)
    }

    # Extract metrics from ASCAT profiles
    QC = ascat.metrics(ascat.bc,ascat.output)

    # Write out segmented regions (including regions with one copy of each allele)
    write.table(ascat.output[["segments"]], file=paste0("${prefix}", ".segments.txt"), sep="\t", quote=F, row.names=F)

    # Write out CNVs in bed format
    cnvs=ascat.output[["segments"]][2:6]
    write.table(cnvs, file=paste0("${prefix}",".cnvs.txt"), sep="\t", quote=F, row.names=F, col.names=T)

    # Write out purity and ploidy info
    summary <- tryCatch({
            matrix(c(ascat.output[["aberrantcellfraction"]], ascat.output[["ploidy"]]), ncol=2, byrow=TRUE)}, error = function(err) {
                # error handler picks up where error was generated
                print(paste("Could not find optimal solution:  ",err))
                return(matrix(c(0,0),nrow=1,ncol=2,byrow = TRUE))
        }
    )
    colnames(summary) <- c("AberrantCellFraction","Ploidy")
    write.table(summary, file=paste0("${prefix}",".purityploidy.txt"), sep="\t", quote=F, row.names=F, col.names=T)

    write.table(QC, file=paste0("${prefix}", ".metrics.txt"), sep="\t", quote=F, row.names=F)

    # Version export
    f <- file("versions.yml","w")
    alleleCounter_version = system(paste("alleleCounter --version"), intern = T)
    ascat_version = as.character(packageVersion('ASCAT'))
    writeLines(paste0('"', "${task.process}", '"', ":"), f)
    writeLines(paste("    ascat:", ascat_version), f)
    writeLines(paste("    alleleCounter:", alleleCounter_version), f)
    close(f)
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.after_correction.gc_rt.test.tumour.germline.png
    touch ${prefix}.after_correction.gc_rt.test.tumour.tumour.png
    touch ${prefix}.before_correction.test.tumour.germline.png
    touch ${prefix}.before_correction.test.tumour.tumour.png
    touch ${prefix}.cnvs.txt
    touch ${prefix}.metrics.txt
    touch ${prefix}.normal_alleleFrequencies_chr21.txt
    touch ${prefix}.normal_alleleFrequencies_chr22.txt
    touch ${prefix}.purityploidy.txt
    touch ${prefix}.segments.txt
    touch ${prefix}.tumour.ASPCF.png
    touch ${prefix}.tumour.sunrise.png
    touch ${prefix}.tumour_alleleFrequencies_chr21.txt
    touch ${prefix}.tumour_alleleFrequencies_chr22.txt
    touch ${prefix}.tumour_normalBAF.txt
    touch ${prefix}.tumour_normalLogR.txt
    touch ${prefix}.tumour_tumourBAF.txt
    touch ${prefix}.tumour_tumourLogR.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-ascat: \$(Rscript -e "library(ASCAT); cat(as.character(packageVersion('ASCAT')))")
        alleleCounter: \$(alleleCounter --version)
    END_VERSIONS
    """
}
