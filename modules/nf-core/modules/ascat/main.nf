process ASCAT {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::ascat=3.0.0 bioconda::cancerit-allelecount-4.3.0": null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-c278c7398beb73294d78639a864352abef2931ce:dfe5aaa885de434adb2b490b68972c5840c6d761-0':
        'quay.io/biocontainers/mulled-v2-c278c7398beb73294d78639a864352abef2931ce:dfe5aaa885de434adb2b490b68972c5840c6d761-0' }"

    input:
    tuple val(meta), path(input_normal), path(index_normal), path(input_tumor), path(index_tumor)
    path(allele_files)
    path(loci_files)

    output:
    tuple val(meta), path("*png"),               emit: png
    tuple val(meta), path("*cnvs.txt"),          emit: cnvs
    tuple val(meta), path("*purityploidy.txt"),  emit: purityploidy
    tuple val(meta), path("*segments.txt"),      emit: segments
    path "versions.yml",                         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args           = task.ext.args        ?: ''
    def prefix         = task.ext.prefix      ?: "${meta.id}"
    def gender         = args.gender          ?  "$args.gender" :        "NULL"
    def genomeVersion  = args.genomeVersion   ?  "$args.genomeVersion" : "NULL"
    def purity         = args.purity          ?  "$args.purity" :        "NULL"
    def ploidy         = args.ploidy          ?  "$args.ploidy" :        "NULL"
    def gc_files       = args.gc_files        ?  "$args.gc_files" :      "NULL"

    def minCounts_arg                    = args.minCounts                     ?  ",minCounts = $args.minCounts" : ""
    def chrom_names_arg                  = args.chrom_names                   ?  ",chrom_names = $args.chrom_names" : ""
    def min_base_qual_arg                = args.min_base_qual                 ?  ",min_base_qual = $args.min_base_qual" : ""
    def min_map_qual_arg                 = args.min_map_qual                  ?  ",min_map_qual = $args.min_map_qual" : ""
    def ref_fasta_arg                    = args.ref_fasta                     ?  ",ref.fasta = '$args.ref_fasta'" : ""
    def skip_allele_counting_tumour_arg  = args.skip_allele_counting_tumour   ?  ",skip_allele_counting_tumour = $args.skip_allele_counting_tumour" : ""
    def skip_allele_counting_normal_arg  = args.skip_allele_counting_normal   ?  ",skip_allele_counting_normal = $args.skip_allele_counting_normal" : ""



    """
    #!/usr/bin/env Rscript
    library(RColorBrewer)
    library(ASCAT)
    options(bitmapType='cairo')


    #prepare from BAM files
    ascat.prepareHTS(
        tumourseqfile = "$input_tumor",
        normalseqfile = "$input_normal",
        tumourname = "Tumour",
        normalname = "Normal",
        allelecounter_exe = "alleleCounter",
        alleles.prefix = "$allele_files",
        loci.prefix = "$loci_files",
        gender = "$gender",
        genomeVersion = "$genomeVersion",
        nthreads = $task.cpus
        $minCounts_arg
        $chrom_names_arg
        $min_base_qual_arg
        $min_map_qual_arg
        $ref_fasta_arg
        $skip_allele_counting_tumour_arg
        $skip_allele_counting_normal_arg
    )


    #Load the data
    ascat.bc = ascat.loadData(
        Tumor_LogR_file = "Tumour_tumourLogR.txt",
        Tumor_BAF_file = "Tumour_normalBAF.txt",
        Germline_LogR_file = "Tumour_normalLogR.txt",
        Germline_BAF_file = "Tumour_normalBAF.txt",
        genomeVersion = "$genomeVersion",
        gender = "$gender"
    )

    #optional GC wave correction
    if(!is.null($gc_files)){
        ascat.bc = ascat.GCcorrect(ascat.bc, $gc_files)
    }

    #Plot the raw data
    ascat.plotRawData(ascat.bc)

    #Segment the data
    ascat.bc = ascat.aspcf(ascat.bc)

    #Plot the segmented data
    ascat.plotSegmentedData(ascat.bc)

    #Run ASCAT to fit every tumor to a model, inferring ploidy, normal cell contamination, and discrete copy numbers
    #If psi and rho are manually set:
    if (!is.null($purity) && !is.null($ploidy)){
        ascat.output <- ascat.runAscat(ascat.bc, gamma=1, rho_manual=$purity, psi_manual=$ploidy)
    } else if(!is.null($purity) && is.null($ploidy)){
        ascat.output <- ascat.runAscat(ascat.bc, gamma=1, rho_manual=$purity)
    } else if(!is.null($ploidy) && is.null($purity)){
        ascat.output <- ascat.runAscat(ascat.bc, gamma=1, psi_manual=$ploidy)
    } else {
        ascat.output <- ascat.runAscat(ascat.bc, gamma=1)
    }

    #Write out segmented regions (including regions with one copy of each allele)
    write.table(ascat.output[["segments"]], file=paste0("$prefix", ".segments.txt"), sep="\t", quote=F, row.names=F)

    #Write out CNVs in bed format
    cnvs=ascat.output[["segments"]][2:6]
    write.table(cnvs, file=paste0("$prefix",".cnvs.txt"), sep="\t", quote=F, row.names=F, col.names=T)

    #Write out purity and ploidy info
    summary <- tryCatch({
            matrix(c(ascat.output[["aberrantcellfraction"]], ascat.output[["ploidy"]]), ncol=2, byrow=TRUE)}, error = function(err) {
                # error handler picks up where error was generated
                print(paste("Could not find optimal solution:  ",err))
                return(matrix(c(0,0),nrow=1,ncol=2,byrow = TRUE))
        }
    )
    colnames(summary) <- c("AberrantCellFraction","Ploidy")
    write.table(summary, file=paste0("$prefix",".purityploidy.txt"), sep="\t", quote=F, row.names=F, col.names=T)

    #version export. Have to hardcode process name and software name because
    #won't run inside an R-block
    version_file_path="versions.yml"
    f <- file(version_file_path,"w")
    writeLines("ASCAT:", f)
    writeLines(" ascat: 3.0.0",f)
    close(f)
    """


    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo stub > ${prefix}.cnvs.txt
    echo stub > ${prefix}.purityploidy.txt
    echo stub > ${prefix}.segments.txt
    echo stub > Tumour.ASCATprofile.png
    echo stub > Tumour.ASPCF.png
    echo stub > Tumour.germline.png
    echo stub > Tumour.rawprofile.png
    echo stub > Tumour.sunrise.png
    echo stub > Tumour.tumour.png

    echo 'ASCAT:' > versions.yml
    echo ' ascat: 3.0.0' >> versions.yml
    """


}
