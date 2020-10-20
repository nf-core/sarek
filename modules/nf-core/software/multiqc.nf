include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

environment = params.enable_conda ? "bioconda::multiqc=1.9" : null
container = "quay.io/biocontainers/multiqc:1.9--py_1"
if (workflow.containerEngine == 'singularity' && !params.pull_docker_container) container = "https://depot.galaxyproject.org/singularity/multiqc:1.9--py_1"

// Has the run name been specified by the user?
// this has the bonus effect of catching both -name and --name
def custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
}

process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode: params.publish_dir_mode

    conda environment
    container container

    input:
        // path software_versions
        path multiqc_config
        path multiqc_custom_config
        val workflow_summary
        path qc_reports

    output:
        path "*multiqc_report.html"
        path "*_data"
        path "multiqc_plots"

    script:
    title = custom_runName ? "--title \"${custom_runName}\"" : ''
    filename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    custom_config_file = params.multiqc_config ? "--config ${multiqc_custom_config}" : ''
    """
    echo '${workflow_summary}' > workflow_summary_mqc.yaml
    multiqc -f ${title} ${filename} ${custom_config_file} .
    """
}
