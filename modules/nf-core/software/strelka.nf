// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

process STRELKA {
    tag "$meta.id"
    
    label 'cpus_max'
    label 'memory_max'
    
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        // TODO nf-core: If a meta map of sample information is NOT provided in "input:" section
        //               change "publish_id:meta.id" to initialise an empty string e.g. "publish_id:''".
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    // TODO nf-core: Fetch "docker pull" address for latest Biocontainer image of software: e.g. https://biocontainers.pro/#/tools/samtools
    //               If required, multi-tool containers may also be available and are usually named to start with "mulled".
    container "quay.io/biocontainers/strelka:2.9.10--0"

    // TODO nf-core: List required Conda packages.
    //               Software MUST be pinned to channel (i.e. "bioconda") and version (i.e. "1.10") as in the example below.
    //               Pinning the build too e.g. "bioconda::samtools=1.10=h9402c20_2" is not currently a requirement.
    conda (params.conda ? "bioconda::strelka=2.9.10" : null)

    input:
    // TODO nf-core: Where applicable all sample-specific information e.g. "id", "single_end", "read_group"
    //               MUST be provided as an input via a Groovy Map called "meta".
    //               This information may not be required in some instances e.g. indexing reference genome files:
    //               https://github.com/nf-core/modules/blob/master/software/bwa/index/main.nf
    // TODO nf-core: Where applicable please provide/convert compressed files as input/output
    //               e.g. "*.fastq.gz" and NOT "*.fastq", "*.bam" and NOT "*.sam" etc.
    tuple val(meta), path(reads)
    // TODO nf-core: List additional required input channels/values here
    val options

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    // TODO nf-core: If meta is provided in "input:" section then it MUST be added to ALL output channels (except version)
    tuple val(meta), path("*.bam"), emit: bam
    // TODO nf-core: List additional required output channels/values here
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def ioptions = initOptions(options)
    // TODO nf-core: If a meta map of sample information is NOT provided in "input:" section delete the line below
    def prefix   = ioptions.suffix ? "${meta.id}${ioptions.suffix}" : "${meta.id}"
    // TODO nf-core: Where possible, a command MUST be provided to obtain the version number of the software e.g. 1.10
    //               If the software is unable to output a version number on the command-line then it can be manually specified
    //               e.g. https://github.com/nf-core/modules/blob/master/software/homer/annotatepeaks/main.nf
    // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "$ioptions.args" variable
    // TODO nf-core: If the tool supports multi-threading then you MUST provide the appropriate parameter
    //               using the Nextflow "task" variable e.g. "--threads $task.cpus"
    // TODO nf-core: Please indent the command appropriately (4 spaces!!) to help with readability ;)
    """
    software tool \\
        $ioptions.args \\
        --threads $task.cpus \\
        $reads \\
        > ${prefix}.bam
    echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' > ${software}.version.txt
    """
}