// Create consensus VCF from bcftools isec sites.txt output
// Simpler alternative to merging numbered VCFs - only captures presence/absence per caller

process CONSENSUS_FROM_SITES {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/47/474a5ea8dc03366b04df884d89aeacc4f8e6d1ad92266888e7a8e7958d07cde8/data'
        : 'community.wave.seqera.io/library/bcftools_htslib:0a3fa2654b52006f'}"

    input:
    tuple val(meta), path(isec_dir)
    // Expects meta.callers = ['caller1', 'caller2', ...] in same order as isec input

    output:
    tuple val(meta), path("${prefix}.vcf.gz")    , emit: vcf
    tuple val(meta), path("${prefix}.vcf.gz.tbi"), emit: tbi
    tuple val("${task.process}"), val('gawk'), eval("awk --version | head -n1 | sed 's/GNU Awk //; s/, .*//'"), emit: versions_gawk, topic: versions
    tuple val("${task.process}"), val('htslib'), eval("tabix --version | head -n1 | sed 's/tabix (htslib) //'"), emit: versions_htslib, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def callers = meta.callers.join(',')
    """
    awk -v callers="${callers}" '
    BEGIN {
        OFS="\\t"
        n=split(callers, c, ",")
        print "##fileformat=VCFv4.2"
        print "##INFO=<ID=CALLERS,Number=.,Type=String,Description=\\"Variant callers that found this variant\\">"
        print "##INFO=<ID=NCALLERS,Number=1,Type=Integer,Description=\\"Number of callers\\">"
        print "#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO"
    }
    {
        split(\$5, bits, "")
        caller_list=""; count=0
        for (i=1; i<=n; i++) {
            if (bits[i]=="1") {
                caller_list = caller_list (caller_list?",":"") c[i]
                count++
            }
        }
        print \$1, \$2, ".", \$3, \$4, ".", ".", "CALLERS="caller_list";NCALLERS="count
    }' ${isec_dir}/sites.txt | bgzip -c > ${prefix}.vcf.gz

    tabix -p vcf ${prefix}.vcf.gz
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    """
}
