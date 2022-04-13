process QUALIMAP_BAMQC {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::qualimap=2.2.2d" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/qualimap:2.2.2d--1' :
        'quay.io/biocontainers/qualimap:2.2.2d--1' }"

    input:
    tuple val(meta), path(bam)
    path gff

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
    qualimap \\
        --java-mem-size=$memory \\
        bamqc \\
        $args \\
        -bam $bam \\
        $regions \\
        -p $strandedness \\
        $collect_pairs \\
        -outdir $prefix \\
        -nt $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qualimap: \$(echo \$(qualimap 2>&1) | sed 's/^.*QualiMap v.//; s/Built.*\$//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    """
    mkdir -p $prefix/css
    mkdir $prefix/images_qualimapReport
    mkdir $prefix/raw_data_qualimapReport
    cd $prefix/css
    touch agogo.css
    touch basic.css
    touch bgtop.png
    touch comment-close.png
    touch doctools.js
    touch down-pressed.png
    touch jquery.js
    touch plus.png
    touch qualimap_logo_small.png
    touch searchtools.js
    touch up.png
    touch websupport.js
    touch ajax-loader.gif
    touch bgfooter.png
    touch comment-bright.png
    touch comment.png
    touch down.png
    touch file.png
    touch minus.png
    touch pygments.css
    touch report.css
    touch underscore.js
    touch up-pressed.png
    cd ../images_qualimapReport/
    touch genome_coverage_0to50_histogram.png
    touch genome_coverage_quotes.png
    touch genome_insert_size_across_reference.png
    touch genome_mapping_quality_histogram.png
    touch genome_uniq_read_starts_histogram.png
    touch genome_coverage_across_reference.png
    touch genome_gc_content_per_window.png
    touch genome_insert_size_histogram.png
    touch genome_reads_clipping_profile.png
    touch genome_coverage_histogram.png
    touch genome_homopolymer_indels.png
    touch genome_mapping_quality_across_reference.png
    touch genome_reads_content_per_read_position.png
    cd ../raw_data_qualimapReport
    touch coverage_across_reference.txt
    touch genome_fraction_coverage.txt
    touch insert_size_histogram.txt
    touch mapped_reads_nucleotide_content.txt
    touch coverage_histogram.txt
    touch homopolymer_indels.txt
    touch mapped_reads_clipping_profile.txt
    touch mapping_quality_across_reference.txt
    touch duplication_rate_histogram.txt
    touch insert_size_across_reference.txt
    touch mapped_reads_gc-content_distribution.txt
    touch mapping_quality_histogram.txt
    cd ../
    touch genome_results.txt
    touch qualimapReport.html
    cd ../

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qualimap: \$(echo \$(qualimap 2>&1) | sed 's/^.*QualiMap v.//; s/Built.*\$//')
    END_VERSIONS
    """
}
