process DEEPVARIANT {
    tag "$meta.id"
    label 'process_medium'


    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used with DeepVariant at the moment. Please use Docker or Singularity containers."
    }

    container "nvcr.io/nvidia/clara/clara-parabricks:4.0.1-1"

    input:
    tuple val(meta), path(input), path(index), path(intervals)
    path(fasta)
    path(fai)

    output:
    tuple val(meta), path("${prefix}.vcf.gz") ,  emit: vcf
    tuple val(meta), path("${prefix}.g.vcf.gz"),  emit: gvcf
    path "versions.yml"               ,  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def regions = intervals ? "--regions ${intervals}" : ""

    """
    /opt/deepvariant/bin/run_deepvariant \\
        --gpus 1 \\
        --ref=${fasta} \\
        --reads=${input} \\
        --output_vcf=${prefix}.vcf.gz \\
        --output_gvcf=${prefix}.g.vcf.gz \\
        ${args} \\
        ${regions} \\
        --num_shards=${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepvariant: \$(echo \$(/opt/deepvariant/bin/run_deepvariant --version) | sed 's/^.*version //; s/ .*\$//' )
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf.gz
    touch ${prefix}.g.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepvariant: \$(echo \$(/opt/deepvariant/bin/run_deepvariant --version) | sed 's/^.*version //; s/ .*\$//' )
    END_VERSIONS
    """
}


"""
docker run -v $(pwd):/workdir -v $(pwd)/results:/outputdir -w /workdir \\
    --gpus all \\
    nvcr.io/nvidia/clara/clara-parabricks:4.0.1-1 \\
    pbrun deepvariant \\
        --ref /workdir/local_refs/GRCh38/GRCh38.primary_assembly.genome_X.fa \\
        --in-bam /workdir/data/ont_fastq/HG002_20220916-1898.converted.cram \\
        --max-read-size-512 \\
        --run-partition \\
        --mode pacbio \\
        --alt-aligned-pileup none \\
        --vsc-min-fraction-indels 0.06 \\
        --channel-insert-size \\
        --channel-gc-content \\
        --pb-model-file /workdir/deepvariant_extra_model_files_v4.0.1-1/A10/deepvariant.eng \\
        --out-variants /outputdir/HG002_20220916.vcf
"""