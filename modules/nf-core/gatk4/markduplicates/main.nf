process GATK4_MARKDUPLICATES {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e3/e3d753d93f57969fe76b8628a8dfcd23ef44bccd08c4ced7089c1f94bf47c89f/data'
        : 'community.wave.seqera.io/library/gatk4_gcnvkernel_htslib_samtools:d3becb6465454c35'}"

    input:
    tuple val(meta), path(bam)
    path fasta
    path fasta_fai

    output:
    tuple val(meta), path("*cram"), emit: cram, optional: true
    tuple val(meta), path("*bam"), emit: bam, optional: true
    tuple val(meta), path("*.crai"), emit: crai, optional: true
    tuple val(meta), path("*.bai"), emit: bai, optional: true
    tuple val(meta), path("*.metrics"), emit: metrics
    tuple val("${task.process}"), val('gatk4'), eval("gatk --version | sed -n '/GATK.*v/s/.*v//p'"), topic: versions, emit: versions_gatk4
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.bam"

    // If the extension is CRAM, then change it to BAM
    prefix_bam = prefix.tokenize('.')[-1] == 'cram' ? "${prefix.substring(0, prefix.lastIndexOf('.'))}.bam" : prefix

    def input_list = bam.collect { bam_ -> "--INPUT ${bam_}" }.join(' ')
    def reference = fasta ? "--REFERENCE_SEQUENCE ${fasta}" : ""

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK MarkDuplicates] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }

    // Using samtools and not Markduplicates to compress to CRAM speeds up computation:
    // https://medium.com/@acarroll.dna/looking-at-trade-offs-in-compression-levels-for-genomics-tools-eec2834e8b94
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        MarkDuplicates \\
        ${input_list} \\
        --OUTPUT ${prefix_bam} \\
        --METRICS_FILE ${prefix}.metrics \\
        --TMP_DIR . \\
        ${reference} \\
        ${args}

    # If cram files are wished as output, the run samtools for conversion
    if [[ ${prefix} == *.cram ]]; then
        samtools view -Ch -T ${fasta} -o ${prefix} ${prefix_bam}
        rm ${prefix_bam}
        samtools index ${prefix}
    fi
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.bam"
    prefix_no_suffix = task.ext.prefix ? prefix.tokenize('.')[0] : "${meta.id}"
    """
    touch ${prefix_no_suffix}.bam
    touch ${prefix_no_suffix}.cram
    touch ${prefix_no_suffix}.cram.crai
    touch ${prefix_no_suffix}.bai
    touch ${prefix}.metrics
    """
}
