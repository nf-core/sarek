process HTSLIB_BGZIPTABIX {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/86/863ca0dbbba30c8367fa4fbd3fa3a84393532fb7b300a5c5c2e70f0dfc475bbf/data'
        : 'community.wave.seqera.io/library/htslib_xz:32f2772a564b3cd2'}"

    input:
    tuple val(meta), path(infile), path(infile_tbi), path(regions)
    val action
    val make_index
    val out_ext

    output:
    tuple val(meta), path("${outfile}"), emit: output
    tuple val(meta), path("${outfile}.{tbi,csi}"), emit: index, optional: true
    // all htslib tools have the same version, we use bgzip
    tuple val("${task.process}"), val('htslib'), eval("bgzip --version | sed '1! d; s/bgzip (htslib) //'"), topic: versions, emit: versions_htslib
    tuple val("${task.process}"), val('xz'), eval("xz --version | sed '1! d; s/xz (XZ Utils) //'"), topic: versions, emit: versions_xz

    when:
    task.ext.when == null || task.ext.when

    script:
    def allowed_actions = ["compress", "decompress"]
    if (action !in allowed_actions) {
        error("htslib/bgziptabix: Invalid action: ${action}. Allowed actions are: ${allowed_actions.join(', ')}")
    }

    if (action == "decompress" && make_index) {
        log.warn("htslib/bgziptabix: Cannot create index when decompressing. Ignoring make_index option.")
    }

    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    outfile = action == "compress" ? (out_ext ? "${prefix}.${out_ext}.gz" : "${prefix}.gz") : (out_ext ? "${prefix}.${out_ext}" : "${prefix}")

    def compress_cmd = action == "compress" ? "bgzip -c ${args} -@ ${task.cpus}" : "cat"
    def bgzip_cmd = action == "compress" ? "[ '\$(basename ${infile})' != '\$(basename ${outfile})' ] && ln -s ${infile} ${outfile}" : "bgzip -c -d ${args} -@ ${task.cpus} ${infile} > ${outfile}"

    def regions_arg = regions ? "-R ${regions}" : ""
    def tabix_cmd = (make_index && !infile_tbi) ? "tabix -@ ${task.cpus} ${regions_arg} ${args2} -f ${outfile}" : ""
    def link_tabix_cmd = make_index && infile_tbi ? "ln -s ${infile_tbi} ${outfile}.${infile_tbi.extension}" : ""
    def uncompressed_cmd = action == "compress" ? "${compress_cmd} ${infile} > ${outfile}" : (infile.getName() == outfile ? "" : "ln -s ${infile} ${outfile}")
    """
    ${link_tabix_cmd}

    FILE_TYPE=\$(htsfile ${infile})

    case "\$FILE_TYPE" in
        *BGZF-compressed*)
            ${bgzip_cmd} ;;
        *gzip-compressed*)
            [ "\$(basename ${infile})" == "\$(basename ${outfile})" ] && echo "Input and output names cannot be the same" && exit 1
            bgzip -d -c -@ ${task.cpus} ${infile} | ${compress_cmd} > ${outfile} ;;
        *bzip2-compressed*)
            bzcat ${infile} | ${compress_cmd} > ${outfile} ;;
        *XZ-compressed*)
            xzcat ${infile} | ${compress_cmd} > ${outfile} ;;
        *)
            ${uncompressed_cmd} ;;
    esac

    ${tabix_cmd}
    """

    stub:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    outfile = action == "compress" ? (out_ext ? "${prefix}.${out_ext}.gz" : "${prefix}.gz") : (out_ext ? "${prefix}.${out_ext}" : "${prefix}")

    def touch_cmd = action == "compress" ? "echo | bgzip -c" : "echo"
    def index_fmt = args2.contains('-C') ? 'csi' : 'tbi'
    def tabix_cmd = make_index ? "touch ${outfile}.${index_fmt}" : ""
    def link_tabix_cmd = make_index && infile_tbi ? "ln -s ${infile_tbi} ${outfile}.${infile_tbi.extension}" : ""
    """
    echo ${args}

    ${touch_cmd} > ${outfile}

    ${tabix_cmd}
    ${link_tabix_cmd}
    """
}
