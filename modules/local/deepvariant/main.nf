// TODO Remove the module declaration
nextflow.enable.dsl=2

// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided as a string i.e. "options.args"
//               where "params.options" is a Groovy Map that MUST be provided via the addParams section of the including workflow.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.

params.options = [:]
options        = initOptions(params.options)

process DEEPVARIANT {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::deepvariant=1.1.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/deepvariant:1.1.0--py36hf3e76ba_2"
    } else {
        // TODO update the bioconda container to work with run_deepvariant.sh script
        // container "quay.io/biocontainers/deepvariant:1.1.0--py36hf3e76ba_2"
        container "google/deepvariant:1.1.0"
    }

    input:
    // TODO nf-core: Where applicable all sample-specific information e.g. "id", "single_end", "read_group"
    //               MUST be provided as an input via a Groovy Map called "meta".
    //               This information may not be required in some instances e.g. indexing reference genome files:
    //               https://github.com/nf-core/modules/blob/master/software/bwa/index/main.nf
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fai

    output:
    tuple val(meta), path("*.vcf*"),  emit: vcf
    tuple val(meta), path("*.g.vcf*"), emit: gvcf
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    /opt/deepvariant/bin/run_deepvariant \
        --ref=${fasta} \
        --reads=${bam} \
        --output_vcf=${meta.id}.vcf.gz \
        --output_gvcf=${meta.id}.g.vcf.gz \
        ${options.args}

    echo \$(/opt/deepvariant/bin/run_deepvariant --version)  > ${software}.version.txt
    """


    stub:

    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    echo "/opt/deepvariant/bin/run_deepvariant \
        --ref=${fasta} \
        --reads=${bam} \
        --output_vcf=${meta.id}.vcf.gz \
        --output_gvcf=${meta.id}.g.vcf.gz"

    echo "${options.args}"

    touch ${meta.id}.vcf.gz
    touch ${meta.id}.g.vcf.gz

    echo \$(/opt/deepvariant/bin/run_deepvariant --version)  > ${software}.version.txt
    """

}


workflow test {

    bam_ch = Channel.of([[id: "HG003.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20"],
                        "${launchDir}/data/HG003.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam",
                        "${launchDir}/data/HG003.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam.bai"]
    )

    fasta_ch = Channel.of("${launchDir}/data/GRCh38_no_alt_analysis_set.fasta")

    fai_ch = Channel.of("${launchDir}/data/GRCh38_no_alt_analysis_set.fasta.fai")

    DEEPVARIANT(bam_ch, fasta_ch, fai_ch )


}
