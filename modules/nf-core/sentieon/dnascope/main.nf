process SENTIEON_DNASCOPE {
    tag "$meta.id"
    label 'process_high'
    label 'sentieon'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a6/a64461f38d76bebea8e21441079e76e663e1168b0c59dafee6ee58440ad8c8ac/data' :
        'community.wave.seqera.io/library/sentieon:202308.03--59589f002351c221' }"

    input:
    tuple val(meta), path(bam), path(bai), path(intervals)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(dbsnp)
    tuple val(meta5), path(dbsnp_tbi)
    tuple val(meta6), path(ml_model)
    val(pcr_indel_model)
    val(emit_vcf)
    val(emit_gvcf)

    output:
    tuple val(meta), path("*.unfiltered.vcf.gz")    , optional:true, emit: vcf   // added the substring ".unfiltered" in the filename of the vcf-files since without that the g.vcf.gz-files were ending up in the vcf-channel
    tuple val(meta), path("*.unfiltered.vcf.gz.tbi"), optional:true, emit: vcf_tbi
    tuple val(meta), path("*.g.vcf.gz")             , optional:true, emit: gvcf   // these output-files have to have the extension ".vcf.gz", otherwise the subsequent GATK-MergeVCFs will fail.
    tuple val(meta), path("*.g.vcf.gz.tbi")         , optional:true, emit: gvcf_tbi
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args                      = task.ext.args                      ?: ''  // options for the driver
    def args2                     = task.ext.args2                     ?: ''  // options for the vcf generation
    def args3                     = task.ext.args3                     ?: ''  // options for the gvcf generation
    def interval                  = intervals                          ? "--interval ${intervals}"               : ''
    def dbsnp_cmd                 = dbsnp                              ? "-d ${dbsnp}"                           : ''
    def model_cmd                 = ml_model                           ? " --model ${ml_model}"                  : ''
    def pcr_indel_model_cmd       = pcr_indel_model                    ? " --pcr_indel_model ${pcr_indel_model}" : ''
    def prefix                    = task.ext.prefix                    ?: "${meta.id}"
    def vcf_cmd                   = ""
    def gvcf_cmd                  = ""
    def base_cmd                  = '--algo DNAscope ' + dbsnp_cmd + ' '

    if (emit_vcf) {  // emit_vcf can be the empty string, 'variant', 'confident' or 'all' but NOT 'gvcf'
        vcf_cmd = base_cmd + args2 + ' ' + model_cmd + pcr_indel_model_cmd + ' --emit_mode ' + emit_vcf + ' ' + prefix + '.unfiltered.vcf.gz'
    }

    if (emit_gvcf) { // emit_gvcf can be either true or false
        gvcf_cmd = base_cmd + args3 + ' ' + model_cmd + pcr_indel_model_cmd + ' --emit_mode gvcf ' + prefix + '.g.vcf.gz'
    }

    def sentieonLicense = secrets.SENTIEON_LICENSE_BASE64 ?
        "export SENTIEON_LICENSE=\$(mktemp);echo -e \"${secrets.SENTIEON_LICENSE_BASE64}\" | base64 -d > \$SENTIEON_LICENSE; " :
        ""
    """
    $sentieonLicense

    sentieon driver $args -r $fasta -t $task.cpus -i $bam $interval $vcf_cmd $gvcf_cmd

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.unfiltered.vcf.gz
    touch ${prefix}.unfiltered.vcf.gz.tbi
    touch ${prefix}.g.vcf.gz
    touch ${prefix}.g.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g" )
    END_VERSIONS
    """
}
