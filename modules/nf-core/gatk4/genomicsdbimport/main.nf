process GATK4_GENOMICSDBIMPORT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.5.0.0--py36hdfd78af_0':
        'biocontainers/gatk4:4.5.0.0--py36hdfd78af_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi), path(interval_file), val(interval_value), path(wspace)
    val   run_intlist
    val   run_updatewspace
    val   input_map

    output:
    tuple val(meta), path("$prefix")        , optional:true, emit: genomicsdb
    tuple val(meta), path("$updated_db")    , optional:true, emit: updatedb
    tuple val(meta), path("*.interval_list"), optional:true, emit: intervallist
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"

    // settings for running default create gendb mode
    input_command = input_map ? "--sample-name-map ${vcf[0]}" : vcf.collect(){"--variant $it"}.join(' ')

    genomicsdb_command = "--genomicsdb-workspace-path ${prefix}"
    interval_command = interval_file ? "--intervals ${interval_file}" : "--intervals ${interval_value}"
    updated_db = ""

    // settings changed for running get intervals list mode if run_intlist is true
    if (run_intlist) {
        genomicsdb_command = "--genomicsdb-update-workspace-path ${wspace}"
        interval_command = "--output-interval-list-to-file ${prefix}.interval_list"
    }

    // settings changed for running update gendb mode. input_command same as default, update_db forces module to emit the updated gendb
    if (run_updatewspace) {
        genomicsdb_command = "--genomicsdb-update-workspace-path ${wspace}"
        interval_command = ''
        updated_db = "${wspace}"
    }

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK GenomicsDBImport] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        GenomicsDBImport \\
        $input_command \\
        $genomicsdb_command \\
        $interval_command \\
        --tmp-dir . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    prefix   = task.ext.prefix ?: "${meta.id}"

    genomicsdb_command = "--genomicsdb-workspace-path ${prefix}"
    interval_command = interval_file ? "--intervals ${interval_file}" : "--intervals ${interval_value}"
    updated_db = ""

    // settings changed for running get intervals list mode if run_intlist is true
    if (run_intlist) {
        genomicsdb_command = "--genomicsdb-update-workspace-path ${wspace}"
        interval_command = "--output-interval-list-to-file ${prefix}.interval_list"
    }

    // settings changed for running update gendb mode. input_command same as default, update_db forces module to emit the updated gendb
    if (run_updatewspace) {
        genomicsdb_command = "--genomicsdb-update-workspace-path ${wspace}"
        interval_command = ''
        updated_db = "${wspace}"
    }

    def stub_genomicsdb = genomicsdb_command == "--genomicsdb-workspace-path ${prefix}" ? "touch ${prefix}" : ""
    def stub_interval   = interval_command == "--output-interval-list-to-file ${prefix}.interval_list" ? "touch ${prefix}.interval_list" : ""
    def stub_update     = updated_db != "" ? "touch ${wspace}" : ""

    """
    ${stub_genomicsdb}
    ${stub_interval}
    ${stub_update}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
