process XENGSORT_INDEX {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/xengsort:2.1.0--pyhdfd78af_1':
        'community.wave.seqera.io/library/htslib_pip_xengsort:39f2d076bf655602'}"

    input:
    path(graft_ref), stageAs: 'graft_ref/*' // allows multiple ref
    path(host_ref) , stageAs: 'host_ref/*'  // allows multiple ref

    output:
    path "xengsort_index", emit: index
    tuple val("${task.process}"), val('xengsort'), eval("xengsort --version"), topic: versions, emit: versions_xengsort_index

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // TODO: Add a separate module to automatically estimate the optimal --nobjects value
    //       for xengsort index based on the size of the reference FASTA files, such as ntCard or KmerEstimate,
    //       as suggested in the xengsort documentation: https://gitlab.com/genomeinformatics/xengsort/-/blob/1ae75047bc76e0daf91c3223245991a8ac39bcfe/README.md
    def nobjects = task.ext.nobjects ?: ''
    def kmersize = task.ext.kmersize ?: ''
    def mask = task.ext.mask ?: ''

    if (!nobjects) {
        error "ERROR: Missing required argument 'nobjects' for xengsort index. Details: https://gitlab.com/genomeinformatics/xengsort/-/blob/1ae75047bc76e0daf91c3223245991a8ac39bcfe/README.md"
    }

    if (kmersize && mask) {
        error "ERROR: Mutually exclusive arguments. Specify only 'kmersize' OR 'mask', not both."
    } else if (!kmersize && !mask) {
        error "ERROR: Missing required argument. Specify either 'kmersize' or 'mask'."
    }

    def kmer_or_mask = kmersize ? "--kmersize ${kmersize}" : "--mask ${mask}"

    """
    mkdir -p xengsort_index

    xengsort index \\
        --index xengsort_index/index \\
        --graft $graft_ref \\
        --host $host_ref \\
        $kmer_or_mask \\
        --nobjects $nobjects \\
        --threads-read ${task.cpus} \\
        --threads-split ${task.cpus} \\
        $args
    """

    stub:
    def args = task.ext.args ?: ''

    """
    echo $args

    mkdir -p xengsort_index
    touch xengsort_index/index.info
    touch xengsort_index/index.hash
    """
}
