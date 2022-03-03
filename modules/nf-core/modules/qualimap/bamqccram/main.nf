process QUALIMAP_BAMQCCRAM {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::qualimap=2.2.2d bioconda::samtools=1.15" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-d3934ca6bb4e61334891ffa2e9a4c87a530e3188:9838874d42d4477d5042782ee019cec9854da7d5-0' :
        'quay.io/biocontainers/mulled-v2-d3934ca6bb4e61334891ffa2e9a4c87a530e3188:9838874d42d4477d5042782ee019cec9854da7d5-0' }"

    input:
    tuple val(meta), path(cram), path(crai)
    path  gff
    path  fasta
    path  fasta_fai

    output:
    tuple val(meta), path("${prefix}"), emit: results
    path  "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"

    def collect_pairs = meta.single_end ? '' : '--collect-overlap-pairs'
    def memory     = task.memory.toGiga() + "G"
    def regions = gff ? "--gff $gff" : ''

    def strandedness = 'non-strand-specific'
    if (meta.strandedness == 'forward') {
        strandedness = 'strand-specific-forward'
    } else if (meta.strandedness == 'reverse') {
        strandedness = 'strand-specific-reverse'
    }
    """
    unset DISPLAY
    mkdir tmp
    export _JAVA_OPTIONS=-Djava.io.tmpdir=./tmp

    samtools view -hb -T ${fasta} ${cram} |
    qualimap \\
        --java-mem-size=$memory \\
        bamqc \\
        $args \\
        -bam /dev/stdin \\
        $regions \\
        -p $strandedness \\
        $collect_pairs \\
        -outdir $prefix \\
        -nt $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qualimap: \$(echo \$(qualimap 2>&1) | sed 's/^.*QualiMap v.//; s/Built.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
