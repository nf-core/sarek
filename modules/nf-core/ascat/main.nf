process ASCAT {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::ascat=3.0.0 bioconda::cancerit-allelecount=4.3.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-c278c7398beb73294d78639a864352abef2931ce:dfe5aaa885de434adb2b490b68972c5840c6d761-0':
        'quay.io/biocontainers/mulled-v2-c278c7398beb73294d78639a864352abef2931ce:dfe5aaa885de434adb2b490b68972c5840c6d761-0' }"

    input:
    tuple val(meta), path(input_normal), path(index_normal), path(input_tumor), path(index_tumor)
    path(allele_files)
    path(loci_files)
    path(bed_file)  // optional
    path(fasta)     // optional
    path(gc_file)   // optional
    path(rt_file)   // optional

    output:
    tuple val(meta), path("*alleleFrequencies_chr*.txt"),      emit: allelefreqs
    tuple val(meta), path("*BAF.txt"),                         emit: bafs
    tuple val(meta), path("*cnvs.txt"),                        emit: cnvs
    tuple val(meta), path("*LogR.txt"),                        emit: logrs
    tuple val(meta), path("*metrics.txt"),                     emit: metrics
    tuple val(meta), path("*png"),                             emit: png
    tuple val(meta), path("*purityploidy.txt"),                emit: purityploidy
    tuple val(meta), path("*segments.txt"),                    emit: segments
    path "versions.yml",                                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args           = task.ext.args        ?: ''
    def prefix         = task.ext.prefix      ?: "${meta.id}"
    def gender         = args.gender          ?  "$args.gender" :        "NULL"
    def genomeVersion  = args.genomeVersion   ?  "$args.genomeVersion" : "NULL"
    def purity         = args.purity          ?  "$args.purity" :        "NULL"
    def ploidy         = args.ploidy          ?  "$args.ploidy" :        "NULL"
    def gc_input       = gc_file              ?  "$gc_file" :            "NULL"
    def rt_input       = rt_file              ?  "$rt_file" :            "NULL"

    def minCounts_arg                    = args.minCounts                     ?  ",minCounts = $args.minCounts" : ""
    def bed_file_arg                     = bed_file                           ?  ",BED_file = '$bed_file'": ""
    def chrom_names_arg                  = args.chrom_names                   ?  ",chrom_names = $args.chrom_names" : ""
    def min_base_qual_arg                = args.min_base_qual                 ?  ",min_base_qual = $args.min_base_qual" : ""
    def min_map_qual_arg                 = args.min_map_qual                  ?  ",min_map_qual = $args.min_map_qual" : ""
    def fasta_arg                        = fasta                              ?  ",ref.fasta = '$fasta'" : ""
    def skip_allele_counting_tumour_arg  = args.skip_allele_counting_tumour   ?  ",skip_allele_counting_tumour = $args.skip_allele_counting_tumour" : ""
    def skip_allele_counting_normal_arg  = args.skip_allele_counting_normal   ?  ",skip_allele_counting_normal = $args.skip_allele_counting_normal" : ""

    """
    #!/usr/bin/env Rscript
    library(RColorBrewer)
    library(ASCAT)
    options(bitmapType='cairo')

    #build prefixes: <abspath_to_files/prefix_chr>
    allele_path = normalizePath("$allele_files")
    allele_prefix = paste0(allele_path, "/", "$allele_files", "_chr")

    loci_path = normalizePath("$loci_files")
    loci_prefix = paste0(loci_path, "/", "$loci_files", "_chr")

    #prepare from BAM files
    ascat.prepareHTS(
        tumourseqfile = "$input_tumor",
        normalseqfile = "$input_normal",
        tumourname = paste0("$prefix", ".tumour"),
        normalname = paste0("$prefix", ".normal"),
        allelecounter_exe = "alleleCounter",
        alleles.prefix = allele_prefix,
        loci.prefix = loci_prefix,
        gender = "$gender",
        genomeVersion = "$genomeVersion",
        nthreads = $task.cpus
        $minCounts_arg
        $bed_file_arg
        $chrom_names_arg
        $min_base_qual_arg
        $min_map_qual_arg
        $fasta_arg
        $skip_allele_counting_tumour_arg
        $skip_allele_counting_normal_arg
    )


    #Load the data
    ascat.bc = ascat.loadData(
        Tumor_LogR_file = paste0("$prefix", ".tumour_tumourLogR.txt"),
        Tumor_BAF_file = paste0("$prefix", ".tumour_tumourBAF.txt"),
        Germline_LogR_file = paste0("$prefix", ".tumour_normalLogR.txt"),
        Germline_BAF_file = paste0("$prefix", ".tumour_normalBAF.txt"),
        genomeVersion = "$genomeVersion",
        gender = "$gender"
    )

    #Plot the raw data
    ascat.plotRawData(ascat.bc, img.prefix = paste0("$prefix", ".before_correction."))

    # optional LogRCorrection
    if("$gc_input" != "NULL") {
        gc_input = paste0(normalizePath("$gc_input"), "/", "$gc_input", ".txt")

        if("$rt_input" != "NULL"){
            rt_input = paste0(normalizePath("$rt_input"), "/", "$rt_input", ".txt")
            ascat.bc = ascat.correctLogR(ascat.bc, GCcontentfile = gc_input, replictimingfile = rt_input)
            #Plot raw data after correction
            ascat.plotRawData(ascat.bc, img.prefix = paste0("$prefix", ".after_correction_gc_rt."))
        }
        else {
            ascat.bc = ascat.correctLogR(ascat.bc, GCcontentfile = gc_input, replictimingfile = $rt_input)
            #Plot raw data after correction
            ascat.plotRawData(ascat.bc, img.prefix = paste0("$prefix", ".after_correction_gc."))
        }
    }

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

    #Extract metrics from ASCAT profiles
    QC = ascat.metrics(ascat.bc,ascat.output)

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

    write.table(QC, file=paste0("$prefix", ".metrics.txt"), sep="\t", quote=F, row.names=F)

    # version export
    f <- file("versions.yml","w")
    alleleCounter_version = system(paste("alleleCounter --version"), intern = T)
    ascat_version = sessionInfo()\$otherPkgs\$ASCAT\$Version
    writeLines(paste0('"', "$task.process", '"', ":"), f)
    writeLines(paste("    alleleCounter:", alleleCounter_version), f)
    writeLines(paste("    ascat:", ascat_version), f)
    close(f)

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo stub > ${prefix}.after_correction.gc_rt.test.tumour.germline.png
    echo stub > ${prefix}.after_correction.gc_rt.test.tumour.tumour.png
    echo stub > ${prefix}.before_correction.test.tumour.germline.png
    echo stub > ${prefix}.before_correction.test.tumour.tumour.png
    echo stub > ${prefix}.cnvs.txt
    echo stub > ${prefix}.metrics.txt
    echo stub > ${prefix}.normal_alleleFrequencies_chr21.txt
    echo stub > ${prefix}.normal_alleleFrequencies_chr22.txt
    echo stub > ${prefix}.purityploidy.txt
    echo stub > ${prefix}.segments.txt
    echo stub > ${prefix}.tumour.ASPCF.png
    echo stub > ${prefix}.tumour.sunrise.png
    echo stub > ${prefix}.tumour_alleleFrequencies_chr21.txt
    echo stub > ${prefix}.tumour_alleleFrequencies_chr22.txt
    echo stub > ${prefix}.tumour_normalBAF.txt
    echo stub > ${prefix}.tumour_normalLogR.txt
    echo stub > ${prefix}.tumour_tumourBAF.txt
    echo stub > ${prefix}.tumour_tumourLogR.txt

    echo "${task.process}:" > versions.yml
    echo ' alleleCounter: 4.3.0' >> versions.yml
    echo ' ascat: 3.0.0' >> versions.yml

    """


}
